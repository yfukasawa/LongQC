#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <zlib.h>
#include <math.h>
#include <inttypes.h>
#include "minimap.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

//YF added.
#include <argp.h>
#include "kalloc.h"
#include "khash.h"
#include "mmpriv.h"
#include "lqutils.h"
#include "kthread.h"
#include "kvec.h"
#include "minimap2-coverage.h"
KHASH_MAP_INIT_INT64(64, uint64_t)

static inline int64_t mm_parse_num(const char *str)
{
    double x;
    char *p;
    x = strtod(str, &p);
    if (*p == 'G' || *p == 'g') x *= 1e9;
    else if (*p == 'M' || *p == 'm') x *= 1e6;
    else if (*p == 'K' || *p == 'k') x *= 1e3;
    return (int64_t)(x + .499);
}

const char *argp_program_version = "minimap2-coverage 0.3 (forked from 2.6-r639-dirty)";
const char *argp_program_bug_address = "<yoshinori.fukasawa@kaust.edu.sa>";
//static char usage[] = "minimap2-coverage [options] <target.seqs> <query.seqs>";

struct arguments
{
    short h_flag, ava_flag, avs_flag, filter_flag, minimizer_cnt_flag;
    int load_n_threads;
    int map_n_threads;
    int min_coverage;
    int num_subsetseq;
    int flag_mm_ava;
    int flag_lq_ava;
    int max_gap;
    int min_cnt;
    int min_score;
    int min_score_med;
    int min_score_good;
    int chain_skip;
    int max_ohang;
    int min_ovlp;
    //int min_match;
    int param_k;
    int param_w;
    double min_ratio;
    uint64_t batch_size;
    char *args[2];
    char *dumpf;
};

static int
parse_opt (int key, char *arg,
           struct argp_state *state)
{
    struct arguments *a = (struct arguments*)state->input;
    switch (key)
    {
        case 'H':
            a->h_flag = 1;
            break;
        case 'k':
            a->param_k = atoi(arg);
            break;
        case 'w':
            a->param_w = atoi(arg);
            break;
        case 'I':
            a->batch_size = mm_parse_num(arg);
            break;
        case 'd':
            a->dumpf = arg;
            break;
        case 'g':
            a->max_gap = atoi(arg);
            break;
        case 'n':
            a->min_cnt = atoi(arg);
            break;
        case 'm':
            a->min_score = atoi(arg);
            break;
        case 'p':
            a->min_score_med = atoi(arg);
            break;
        case 'q':
            a->min_score_good = atoi(arg);
            break;
        case 's':
            a->chain_skip = atoi(arg);
            break;
        case 'X':
            a->ava_flag = 1;
            break;
        case 'Y':
            a->avs_flag = 1;
            break;
        case 'a':
            a->max_ohang = atoi(arg);
            break;
        case 'l':
            a->min_ovlp = atoi(arg);
            break;
        case 'c':
            a->min_coverage = atoi(arg);
            break;
        case 'r':
            a->min_ratio = atof(arg);
            break;
        case 'f':
            a->filter_flag = 1;
            break;
        case 'z':
            a->minimizer_cnt_flag = 1;
            break;
        //case 'e':
        //    a->min_match = atoi(arg);
        //    break;
        case 'u':
            a->num_subsetseq = atoi(arg);
            //printf ("%.4f\n", a->min_ratio);
            break;
        case 't':
            a->map_n_threads = atoi(arg);
            break;
        case ARGP_KEY_ARG:
            if (state->arg_num >= 2)
            /* Too many arguments. */
                argp_usage (state);
            a->args[state->arg_num] = arg;
            break;
        case ARGP_KEY_INIT:
            a->param_k = a->param_w = 0;
            a->ava_flag = a->avs_flag = a->h_flag = a->filter_flag = 0;
            a->load_n_threads = 3;
            a->map_n_threads  = 1;
            a->min_coverage = a->num_subsetseq = -1;
            a->max_gap = a->min_cnt = a->min_score = a->min_score_med = a->min_score_good = 0;
            a->chain_skip = a->max_ohang = a->min_ovlp = -1;
            a->min_ratio  = 0.0; //a->min_match = 3;
            a->dumpf = 0; a->args[0] = 0; a->args[1] = 0;
            a->batch_size = 0;
            break;
        case ARGP_KEY_END:
            if (a->dumpf == 0 && state->arg_num < 2)
            /* Not enough arguments. */
                argp_usage (state);
            break;
        default:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

struct argp_option options[] =
{
    { 0, 0, 0, 0, "Indexing options:", 1},
    { "homopolymer",       'H',     0,  0, "use homopolymer-compressed k-mer"},
    { "k-mer",             'k', "INT",  0, "k-mer size (no larger than 28)"},
    { "window",            'w', "INT",  0, "minizer window size"},
    { "index-size",        'I', "STRING",  0, "split index for every ~NUM input bases"},
    { "dump-index",        'd', "FILE", 0, "dump index to FILE"},
    { 0, 0, 0, 0, "Mapping options", 2},
    { "max-gap-length",    'g', "INT",  0, "stop chain enlongation if there are no minimizers in INT-bp"},
    { "min-cnt",           'n', "INT",  0, "minimal number of minimizers on a chain"},
    { "min-score",         'm', "INT",  0, "minimal chaining score (matching bases minus log gap penalty)"},
    { "min-score-t2",      'p', "INT",  0, "chaining score to filter bad reads and must be >= min-score"},
    { "min-score-t3",      'q', "INT",  0, "chaining score to filter bad and some spurious reads/overlap and must be >= min-score-t2"},
    { "max-chain-skip",    's', "INT",  OPTION_HIDDEN, "behaves in the same way of that in minimap2"},
    { "skip-self-ava",     'X',     0,  0, "skip self and dual mappings (for the all-vs-all mode)"},
    { "skip-self",         'Y',     0,  0, "skip self mappings (for the all-vs-subfraction mode)"},
    { 0, 0, 0, 0, "Filtering options", 3},
    { "max-overhang",      'a', "INT",  0, "stop chain enlongation if there are no minimizers in INT-bp"},
    { "min-overlap-len",   'l', "INT",  0, "minimal number of minimizers on a chain"},
    { "min-coverage",      'c', "INT",  0, "minimal chaining score (matching bases minus log gap penalty)"},
    { "min-overlap-ratio", 'r', "NUM",  0, "minimum coverage ratio. overlap length is compared with overhang length."},
    //{ "min-match-kmer",    'e', "INT",  0, "minimizers having more than this num are considered as error-free k-mer"},
    { 0, 0, 0, 0, "Misc options", 4},
    { "num-subset",        'u', "INT",  0, "number of sequences for query subset"},
    { "threads",           't', "INT",  0, "number of threads"},
    { "minimizer-cnt",     'z',     0,  0, "count minimizer apperance for all queries"},
    { "filter",            'f',     0,  0, "read filtering mode. Can be used to filter out known contaminants (e.g. spiked-in control reads)."},
    { 0 }
};

static struct argp argp = { options, parse_opt, "reference reads", "minimap2-coverage is a forked minimap2 program focusing on coverage, and tuned for LongQC.\v"" " };

// for debug
void printb(unsigned int v) {
    unsigned int mask = (int)1 << (sizeof(v) * CHAR_BIT - 1);
    do putchar(mask & v ? '1' : '0');
    while (mask >>= 1);
}

int main(int argc, char *argv[])
{
    int offset_q, i, j, k, buf;
    buf = 10000;
    
    kseq_t *ks;
    mm128_v mv;
    lq_minimizer_cnt_v m_cnts;
    gzFile f;

    struct arguments a;
    argp_parse (&argp, argc, argv, 0, 0, &a);
    
    mm_idxopt_t iopt;
    mm_mapopt_t mopt;
    lq_fltopt_t fopt;
    
    mm_idxopt_init(&iopt);
    mm_mapopt_init(&mopt);
    
    mm_verbose = 2; // disable message output to stderr
    mm_realtime0 = realtime();
    
    if (a.ava_flag != 0 && a.avs_flag != 0){
        fprintf(stderr, "Error: -X and -Y are mutually exclusive\n");
        return 1;
    } else if (a.ava_flag == 0 && a.avs_flag == 0 && a.dumpf == 0) {
        fprintf(stderr, "Error: Choose either -X (all-vs-all) or -Y (all-vs-sub)\n");
        return 1;
    } else if(a.ava_flag) {
        mopt.flag |= MM_F_NO_SELF; //YF memo: this doesn't work in single thread version. bizzare...
        mopt.flag |= MM_F_AVA;
    } else if(a.avs_flag) {
        mopt.flag |= MM_F_NO_SELF; //YF memo: this doesn't work in single thread version. bizzare...
        mopt.flag |= LQ_F_AVA;
    }
    
    if(a.h_flag == 1){
        iopt.flag |= MM_I_HPC;
    } else if (a.h_flag == 0) {
        fprintf(stderr, "Homopolymer compression is not applied.\n");
    } else {
        fprintf(stderr, "Error: unexpected value was given to -H option. \n");
        return 1;
    }
    
    if(a.param_k == 0){
        fprintf(stderr, "Warning: Apply default k=12 instead. \n");
        iopt.k = 12;
    } else {
        iopt.k = a.param_k;
    }
    
    if(a.param_w == 0){
        fprintf(stderr, "Warning: Apply default w=5 instead. \n");
        iopt.w = 5;
    } else {
        iopt.w = a.param_w;
    }
    
    if(a.batch_size == 0){
        fprintf(stderr, "Warning: Apply default I=4G instead. \n");
    } else {
        iopt.batch_size = a.batch_size;
    }
    
    // open index reader
    mm_idx_reader_t *r = mm_idx_reader_open(a.args[0], &iopt, a.dumpf);
    mm_idx_t *mi;
    
    if (r == 0) {
        fprintf(stderr, "ERROR: failed to open file '%s'\n", a.args[0]);
        return 1;
    }
    
    /*
     if (!r->is_idx && a.dumpf == 0 && a.args[1] == 0) {
     fprintf(stderr, "ERROR: input is missing. please specify a query file or option -d to keep the index\n");
     mm_idx_reader_close(r);
     return 1;
     }
     */

    
    if(a.min_coverage == -1){
        fprintf(stderr, "Warning: Apply default c=3 instead. \n");
        fopt.min_coverage = 3;
    } else {
        fopt.min_coverage = a.min_coverage;
    }
    
    if(a.num_subsetseq == -1){
        fprintf(stderr, "Warning: -s shouldn't be zero. Apply default s=100000 instead.\n");
        a.num_subsetseq = 100000;
    }

    if(a.max_gap == 0){
        fprintf(stderr,"Warning: Apply default g=10000 instead.\n");
        mopt.max_gap = 10000;
    } else {
        mopt.max_gap = a.max_gap;
    }
    
    if(a.min_cnt == 0){
        fprintf(stderr, "Warning: Apply default n=3 instead.\n");
        mopt.min_cnt = 3;
    } else {
        mopt.min_cnt = a.min_cnt;
    }
    
    if(a.min_score == 0){
        fprintf(stderr, "Warning: Apply default m=40 instead.\n");
        mopt.min_chain_score = 40;
    } else {
        mopt.min_chain_score = a.min_score;
    }
    
    // YF memo: because of the default, this has to come after mopt.min_chain_score loading.
    if(a.min_score_med == 0){
        fprintf(stderr, "Warning: Apply default p=m instead.\n");
        a.min_score_med = mopt.min_chain_score;
    }
    
    if(a.min_score_good == 0){
        fprintf(stderr, "Warning: Apply default q=m instead.\n");
        a.min_score_good = mopt.min_chain_score;
    }
    
    if(a.min_score_med != 0 && mopt.min_chain_score != 0){
        if(a.min_score_med < mopt.min_chain_score){
            fprintf(stderr, "Error: -p must be larger than or equal to -m.\n");
            exit(1);
        }
    }
    
    if(a.min_score_good != 0 && a.min_score_med!= 0 && mopt.min_chain_score != 0){
        if(a.min_score_good < mopt.min_chain_score || a.min_score_good < a.min_score_med){
            fprintf(stderr, "Error: -q must be larger than or equal to -m and -p.\n");
            exit(1);
        }
    }
    
    if(a.min_score_good != 0 && a.min_score_med!= 0 && mopt.min_chain_score == 0){
        if(a.min_score_good < a.min_score_med){
            fprintf(stderr, "Error: -q must be larger than or equal to -p.\n");
            exit(1);
        }
    }
    
    if(a.min_score_good != 0 && mopt.min_chain_score != 0){
        if(a.min_score_good < mopt.min_chain_score){
            fprintf(stderr, "Error: -q must be larger than or equal to -m.\n");
            exit(1);
        }
    }
    
    if(a.chain_skip == -1){
        fprintf(stderr, "Warning: Apply default s=25 instead.\n");
        mopt.max_chain_skip  = 25;
    } else {
        mopt.max_chain_skip  = a.chain_skip;
    }
    
    if(a.max_ohang == -1){
        fprintf(stderr, "Warning: Apply default a=2000 instead.\n");
        fopt.max_overhang = 2000;
    } else {
        fopt.max_overhang = a.max_ohang;
    }
    
    if(a.min_ovlp == -1){
        fprintf(stderr, "Warning: Apply default l=1000 instead.\n");
        fopt.min_ovlp = 1000;
    } else {
        fopt.min_ovlp = a.min_ovlp;
    }
    
    if(a.min_ratio == 0.0){
        fprintf(stderr, "Warning: Apply default r=0.4 instead.\n");
        fopt.min_ratio = 0.4;
    } else {
        fopt.min_ratio = a.min_ratio;
    }
    
    mm_verbose = 3;
    
    fprintf(stderr, "=== Parameters are listed below === \n");
    if(a.args[0] && a.args[1] == NULL){
        fprintf(stderr, "Inputs is target: %s\n", a.args[0]);
    } else if(a.args[0] && a.args[1]) {
        fprintf(stderr, "Inputs are target: %s, query: %s\n", a.args[0], a.args[1]);
    } else {
        fprintf(stderr, "Error: weird inputs.\n");
    }
    fprintf(stderr, "kmer %d, window %d, index loading size %"PRIu64"\n", iopt.k, iopt.w, iopt.batch_size);
    fprintf(stderr, "min-score %d, min-score-med %d, min-score-good %d, max-gap %d, min-cnt %d\n", mopt.min_chain_score, a.min_score_med, a.min_score_good, mopt.max_gap, mopt.min_cnt);
    fprintf(stderr, "Homo-polymer compression: %d, Filtering: %d, minimizer-count: %d\n", a.h_flag, a.filter_flag, a.minimizer_cnt_flag);
    fprintf(stderr, "max-overhang %d, min-overlaplen %d, min-overapratio %.2f\n", fopt.max_overhang, fopt.min_ovlp, fopt.min_ratio);
    fprintf(stderr, "num of threads %d, num of query seqs %d\n===\n", a.map_n_threads, a.num_subsetseq);
    
    if(a.args[1] != NULL){
        // init -> independent file
        f = gzopen(a.args[1], "r");
        ks = kseq_init(f);
        assert(f);

        kv_init(mv);
        kv_init(m_cnts);
        
        i = 0;
        void *km = km_init();
        uint64_t tot = 0;
        while (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence
            mm_sketch(km, ks->seq.s, ks->seq.l, iopt.w, iopt.k, i++, a.h_flag, &mv);
            //printf("[DEBUG] The number of minimisers for %s: %zu\n", ks->name.s, mv.n);
            m_array a;
            a.a = (uint16_t*)calloc(mv.n, sizeof(uint16_t));
            a.n = mv.n; tot += mv.n;
            kv_push(m_array, 0, m_cnts, a);
            kfree(km, mv.a);
            kv_init(mv);
        }
        km_destroy(km);
        fprintf(stderr, "[M::%s::%.3f*%.2f] loaded %d sequence(s). Total m_cnt: %llu (Peak RSS: %.3f GB) \n",
                __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), i, tot, peakrss() / 1024.0 / 1024.0 / 1024.0);
    }
    
    // coverage counter
    // YF memo: each thread access different member (read), and afaik two threads won't access the same one simutaneously.
    uint64_t *lambdas;  lambdas  = (uint64_t*) calloc(m_cnts.n, sizeof(uint64_t)); // chain coverage
    uint64_t *lambdas2; lambdas2 = (uint64_t*) calloc(m_cnts.n, sizeof(uint64_t)); // good chain coverage
    float *avg_ks;        avg_ks = (float*) calloc(m_cnts.n, sizeof(float));
    lq_subcoords_v **ovlp_coords;
    ovlp_coords = (lq_subcoords_v**) malloc(m_cnts.n * sizeof(lq_subcoords_v*));
    for(i=0; i<m_cnts.n; ++i){
        ovlp_coords[i] = (lq_subcoords_v *)malloc(sizeof(lq_subcoords_v));
        kv_init(*ovlp_coords[i]);
        kv_resize(lq_subcoords_t, 0, *ovlp_coords[i], buf); //YF memo: make some buffer to avoid memory fragmentation
    }
    
    // index can be read several times due to the memory size specified -I option. see iopt.batch_size
    // That's why this has to be in a while loop.
    
    offset_q = 0;
    while ((mi = mm_idx_reader_read(r, a.load_n_threads)) != 0) {
        mm_mapopt_update(&mopt, mi); // this sets the maximum minimizer occurrence; TODO: set a better default in mm_mapopt_init()!
        fprintf(stderr, "[M::%s::%.3f*%.2f] loaded/built the index for %d target sequence(s)\n",
                __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), mi->n_seq);
        lq_map_file(mi, a.args[1] == NULL ? NULL : a.args[1], &mopt, &fopt, a.min_score_med, a.min_score_good, a.map_n_threads, lambdas, lambdas2,
                    &m_cnts, ovlp_coords, avg_ks, &offset_q);
        offset_q   = 0; // all of queries were visited by every chunnk of db.
        mm_idx_destroy(mi);
    }
    
    if(a.args[1] == NULL){
        free(lambdas);
        free(lambdas2);
        free(ovlp_coords);
        free(avg_ks);
        mm_idx_reader_close(r); // close the index reader
        fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - mm_realtime0, cputime());
        return 0;
    }
    
    // open query file for reading; you may use your favorite FASTA/Q parser
    f = gzopen(a.args[1], "r");
    assert(f);
    ks = kseq_init(f);
    i = 0;
    uint64_t tot_m, tot_m_cnt;
    tot_m = tot_m_cnt = 0; //debug
    
    if(a.minimizer_cnt_flag){
        // vars to get minimizers. used for paper revision
        mm128_v minis;
        int h_ret;
        void *km = km_init();
        khiter_t itr;
        uint64_t *b; int32_t n_b; uint64_t tot_cnt;
        khash_t(64) *h = kh_init(64);
        kh_resize(64, h, 0);
        
        tot_cnt = 0;
        while (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence
            m_array mv = m_cnts.a[i];
            kv_init(minis);
            mm_sketch(km, ks->seq.s, ks->seq.l, iopt.w, iopt.k, i, a.h_flag, &minis); // YF memo: redundant...
            //printf("[DEBUG] assert %d: %zu\n", mv.n, minis.n);
            assert(mv.n == minis.n);
            for(j=0; j<mv.n; j++){
                /*
                 printf("[DEBUG] %"PRIu64"(i=%d, reverse=%d): %d\n", minis.a[j].x>>8,
                 (uint32_t) minis.a[j].y >>1, (uint32_t) minis.a[j].y&1,
                 mv.a[j]);
                 */
                itr = kh_put(64, h, minis.a[j].x>>8, &h_ret);
                if(h_ret == 0){
                    kh_val(h, itr) = kh_val(h, itr) + (uint64_t) mv.a[j];
                    //printf("[DEBUG] hash value: %"PRIu64"\n", kh_val(h, itr));
                } else if (h_ret < 0){
                    fprintf(stderr, "Error: operation failed.\n"); // some error happend.
                } else if (h_ret == 1) {
                    kh_val(h, itr) = (uint64_t) mv.a[j];
                    //printf("[DEBUG] hash value: %"PRIu64"\n", kh_val(h, itr));
                } else {
                    // return code == 2 should not happen.
                    fprintf(stderr, "Error: unexpected deleion happend.\n"); // some error happend.
                }
            }
            // end of loop processing
            kfree(km, minis.a);
            i++;
        }
        
        n_b = 0;
        b = (uint64_t*)kmalloc(km, kh_size(h)*sizeof(uint64_t));
        for (itr = kh_begin(h); itr != kh_end(h); ++itr)
            if (kh_exist(h, itr))
                b[n_b++] = kh_value(h, itr);
        
        // for showing
        radix_sort_64(b, b + n_b);
        for (j = 0; j < n_b>>1; ++j) { // reverse
            uint64_t t = b[j];
            b[j] = b[n_b - j - 1], b[n_b - j - 1] = t;
        }
        
        //for(j = 0; j < n_b; ++j)
        //    printf("[dev] minimizer %d cnt: %"PRIu64"\n", j, b[j]);
        
        kfree(km, b);
        km_destroy(km);
        kseq_destroy(ks);
        i = 0;
        f = gzopen(a.args[1], "r");
        assert(f);
        ks = kseq_init(f);
    }
    
    while (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence
        lq_subcoords_v regs  = {0,0,0};
        lq_subcoords_v mregs = {0,0,0};
        lq_subcoords_v* v = ovlp_coords[i];
        m_array mv = m_cnts.a[i];
        
        tot_m += v->n; tot_m_cnt += mv.n;
        int32_t n_match = 0; double div = 0.0;
        uint32_t sum = 0.0;
        for(j=0; j<mv.n; j++){
            sum += mv.a[j];
        }
        
        sum /= mv.n;
        for(j=0; j<mv.n; j++){
            if(mv.a[j] > sum) n_match++;
        }
        //printf("[DEBUG] ave. %d, nmatch %d, %.3f\n", sum, n_match, logf((float)mv.n / n_match));
        div = n_match > 0 ? logf((float)mv.n / n_match) / avg_ks[i] : 1.0;
        
        compute_reliable_region(v, fopt.min_coverage, &regs, &mregs);
        
        if(regs.n > 0){
            uint32_t tot = 0;
            char *coords; coords = (char *) malloc(sizeof(char) * (66 *v->n)); coords[0] = '\0';
            for(j=0; j<regs.n; j++){
                if(j>0)
                    strcat(coords, ",");
                char c[66] = {'\0'}; sprintf(c, "%d-%d", regs.a[j].start, regs.a[j].end);
                strcat(coords, c);
                tot += regs.a[j].end - regs.a[j].start;
            }
            if (mregs.n > 0){
                char *mcoords; mcoords = (char *) malloc(sizeof(char) * (66 *v->n)); mcoords[0] = '\0';
                for(k=0; k<mregs.n; k++){
                    if(k>0)
                        strcat(mcoords, ",");
                    char mc[66] = {'\0'}; sprintf(mc, "%d-%d", mregs.a[k].start, mregs.a[k].end);
                    strcat(mcoords, mc);
                }
                
                if(a.filter_flag){
                    printf("%s\t%d\t%"PRIu64"\t%s\t%s\t%.3f\t%.3f\t%.3f\t0.0\n", ks->name.s, ks->seq.l, lambdas[i],
                           coords, mcoords, (double) tot/ks->seq.l, meanQ(ks->qual.s, ks->qual.l), div);
                } else {
                    printf("%s\t%d\t%"PRIu64"\t%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\n", ks->name.s, ks->seq.l, lambdas[i],
                           coords, mcoords, (double) lambdas[i]/tot, meanQ(ks->qual.s, ks->qual.l), div, (double) lambdas2[i]/tot);
                }
                
            } else {
                if(a.filter_flag){
                    printf("%s\t%d\t%"PRIu64"\t%s\t0\t%.3f\t%.3f\t%.3f\t0.0\n", ks->name.s, ks->seq.l, lambdas[i],
                           coords, (double) tot/ks->seq.l, meanQ(ks->qual.s, ks->qual.l), div);
                } else {
                    printf("%s\t%d\t%"PRIu64"\t%s\t0\t%.3f\t%.3f\t%.3f\t%.3f\n", ks->name.s, ks->seq.l, lambdas[i],
                           coords, (double) lambdas[i]/tot, meanQ(ks->qual.s, ks->qual.l), div, (double) lambdas2[i]/tot);
                }
            }
        } else {
            printf("%s\t%d\t%"PRIu64"\t0\t0\t0.0\t%.3f\t%.3f\t0.0\n", ks->name.s, ks->seq.l, lambdas[i], meanQ(ks->qual.s, ks->qual.l), div);
        }
        if(regs.n > 0)
            kv_destroy(regs);
        if(mregs.n > 0)
            kv_destroy(mregs);
        if(v->n > 0)
            kv_destroy(*v);
        free(v);
        
        if(mv.n>0)
            free(mv.a);
        i++;
    }
    //printf("[DEBUG] total array size for coords: %"PRIu64", total array size for cnt: %"PRIu64"\n", tot_m, tot_m_cnt);
    kseq_destroy(ks);
    mm_idx_reader_close(r); // close the index reader
    gzclose(f);
    
//    // open query file for reading; you may use your favorite FASTA/Q parser
//    f = gzopen(a.args[1], "r");
//    assert(f);
//    ks = kseq_init(f);
//    i = 0;
//    uint64_t tot_m, tot_m_cnt;
//    tot_m = tot_m_cnt = 0; //debug
//    while (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence
//        lq_subcoords_v regs  = {0,0,0};
//        lq_subcoords_v mregs = {0,0,0};
//        lq_subcoords_v gregs = {0,0,0};
//        lq_subcoords_v* v = ovlp_coords[i];
//        m_array mv = m_cnts.a[i];
//        tot_m += v->m; tot_m_cnt += mv.n;
//        int32_t n_match = 0; double div = 0.0;
//        uint32_t sum = 0.0;
//        //fprintf(stderr, "debug-match-cnt %s: ", ks->name.s);
//        for(j=0; j<mv.n; j++){
//            sum += mv.a[j];
//        }
//        sum /= mv.n;
//        for(j=0; j<mv.n; j++){
//            if(mv.a[j] > sum) n_match++;
//            //fprintf(stderr, "%d, ", mv->a[j]);
//        }
//        //fprintf(stderr, "ave. %d, nmatch %d\n", sum, n_match);
//        div = n_match > 0 ? logf((float)mv.n / n_match) / avg_ks[i] : 1.0;
//        
//        compute_reliable_region(v, fopt.min_coverage, &regs, &mregs, &gregs);
//        
//        if(regs.n > 0){
//            uint32_t tot = 0;
//            char *coords; coords = (char *) malloc(sizeof(char) * (66 *v->n)); coords[0] = '\0';
//            for(j=0; j<regs.n; j++){
//                if(j>0)
//                    strcat(coords, ",");
//                char c[66] = {'\0'}; sprintf(c, "%d-%d", regs.a[j].start, regs.a[j].end);
//                strcat(coords, c);
//                tot += regs.a[j].end - regs.a[j].start;
//            }
//            if (mregs.n > 0){
//                char *mcoords; mcoords = (char *) malloc(sizeof(char) * (66 *v->n)); mcoords[0] = '\0';
//                for(k=0; k<mregs.n; k++){
//                    if(k>0)
//                        strcat(mcoords, ",");
//                    char mc[66] = {'\0'}; sprintf(mc, "%d-%d", mregs.a[k].start, mregs.a[k].end);
//                    strcat(mcoords, mc);
//                }
//                if (gregs.n > 0){
//                    char *gcoords; gcoords = (char *) malloc(sizeof(char) * (66 *v->n)); gcoords[0] = '\0';
//                    for(l=0; l<gregs.n; l++){
//                        if(l>0)
//                            strcat(gcoords, ",");
//                        char gc[66] = {'\0'}; sprintf(gc, "%d-%d", gregs.a[l].start, gregs.a[l].end);
//                        strcat(gcoords, gc);
//                    }
//                    if(a.filter_flag){
//                        printf("%s\t%d\t%"PRIu64"\t%s\t%s\t%s\t%.3f\t%.3f\t%.3f\t0.0\n", ks->name.s, ks->seq.l, lambdas[i],
//                               coords, mcoords, gcoords, (double) tot/ks->seq.l, meanQ(ks->qual.s, ks->qual.l), div);
//                    } else {
//                        printf("%s\t%d\t%"PRIu64"\t%s\t%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\n", ks->name.s, ks->seq.l, lambdas[i],
//                               coords, mcoords, gcoords, (double) lambdas[i]/tot, meanQ(ks->qual.s, ks->qual.l), div, (double) lambdas2[i]/tot);
//                    }
//                } else {
//                    if(a.filter_flag){
//                        printf("%s\t%d\t%"PRIu64"\t%s\t%s\t0\t%.3f\t%.3f\t%.3f\t0.0\n", ks->name.s, ks->seq.l, lambdas[i],
//                               coords, mcoords, (double) tot/ks->seq.l, meanQ(ks->qual.s, ks->qual.l), div);
//                    } else {
//                        printf("%s\t%d\t%"PRIu64"\t%s\t%s\t0\t%.3f\t%.3f\t%.3f\t%.3f\n", ks->name.s, ks->seq.l, lambdas[i],
//                               coords, mcoords, (double) lambdas[i]/tot, meanQ(ks->qual.s, ks->qual.l), div, (double) lambdas2[i]/tot);
//                    }
//                }
//            } else {
//                if(a.filter_flag){
//                    printf("%s\t%d\t%"PRIu64"\t%s\t0\t0\t%.3f\t%.3f\t%.3f\t0.0\n", ks->name.s, ks->seq.l, lambdas[i],
//                           coords, (double) tot/ks->seq.l, meanQ(ks->qual.s, ks->qual.l), div);
//                } else {
//                    printf("%s\t%d\t%"PRIu64"\t%s\t0\t0\t%.3f\t%.3f\t%.3f\t%.3f\n", ks->name.s, ks->seq.l, lambdas[i],
//                           coords, (double) lambdas[i]/tot, meanQ(ks->qual.s, ks->qual.l), div, (double) lambdas2[i]/tot);
//                }
//            }
//        } else {
//            printf("%s\t%d\t%"PRIu64"\t0\t0\t0\t0.0\t%.3f\t%.3f\t0.0\n", ks->name.s, ks->seq.l, lambdas[i], meanQ(ks->qual.s, ks->qual.l), div);
//        }
//        if(regs.n > 0)
//            kv_destroy(regs);
//        if(mregs.n > 0)
//            kv_destroy(mregs);
//        if(gregs.n > 0)
//            kv_destroy(gregs);
//        if(v->n > 0)
//            kv_destroy(*v);
//        free(v);
//        
//        if(mv.n>0)
//            free(mv.a);
//        i++;
//    }
//    printf("[DEBUG] total array size for coords: %"PRIu64", total array size for cnt: %"PRIu64"\n", tot_m, tot_m_cnt);
//    kseq_destroy(ks);
//    mm_idx_reader_close(r); // close the index reader
//    gzclose(f);
    
    free(lambdas);
    free(lambdas2);
    free(m_cnts.a);
    free(ovlp_coords);
    free(avg_ks);

    fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - mm_realtime0, cputime());
    return 0;
}
