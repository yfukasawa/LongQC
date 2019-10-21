//
//  kmer-dump.c
//  mm2_cov_xcode
//
//  Created by Yoshinori Fukasawa on 19/10/2019.
//  Copyright © 2019 Yoshinori Fukasawa. All rights reserved.
//

#include "kmer-dump.h"

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <zlib.h>
//#include <math.h>
#include <inttypes.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#include <argp.h>
#include "kalloc.h"
#include "khash.h"
#include "mmpriv.h"
//#include "lqutils.h"
#include "kthread.h"
#include "kvec.h"
#include "minimap.h"
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


const char *argp_program_version = "kmer-dump 0.1";
const char *argp_program_bug_address = "<yoshinori.fukasawa@kaust.edu.sa>";
//static char usage[] = "minimap2-coverage [options] <target.seqs> <query.seqs>";

struct arguments
{
    short h_flag, ava_flag, avs_flag, filter_flag;
    int load_n_threads;
    int map_n_threads;
    int param_k;
    int param_w;
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
            a->h_flag = 0;
            a->load_n_threads = 3;
            a->map_n_threads  = 1;
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
    { "index-size",        'I', "STRING", 0, "split index for every ~NUM input bases"},
    { "dump-index",        'd', "FILE", 0, "dump index to FILE"},
    { 0, 0, 0, 0, "Misc options", 2},
    { "threads",           't', "INT",  0, "number of threads"},
    { 0 }
};

static struct argp argp = { options, parse_opt, "reference reads", "kmer-dump is a program dumping occurences of k-mer in the dataset.\v"" " };

// for debug
void printb(unsigned int v) {
    unsigned int mask = (int)1 << (sizeof(v) * CHAR_BIT - 1);
    do putchar(mask & v ? '1' : '0');
    while (mask >>= 1);
}

int main(int argc, char *argv[])
{
    int offset_q, i, j;
    kseq_t *ks;
    mm128_v mv;
    lq_minimizer_cnt_v m_cnts;
    gzFile f;

    struct arguments a;
    argp_parse (&argp, argc, argv, 0, 0, &a);
    
    mm_idxopt_t iopt;
    mm_mapopt_t mopt;
    
    mm_idxopt_init(&iopt);
    mm_mapopt_init(&mopt);
    
    mm_verbose = 2; // disable message output to stderr
    mm_realtime0 = realtime();

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
    fprintf(stderr, "Homo-polymer compression: %d\n", a.h_flag);
    fprintf(stderr, "num of threads %d\n===\n", a.map_n_threads);
    
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
        fprintf(stderr, "[M::%s::%.3f*%.2f] loaded %d sequence(s). Total m_cnt: %"PRIu64" (Peak RSS: %.3f GB) \n",
                __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), i, tot, peakrss() / 1024.0 / 1024.0 / 1024.0);
    }

    // index can be read several times due to the memory size specified -I option. see iopt.batch_size
    // That's why this has to be in a while loop.
    offset_q = 0;
    while ((mi = mm_idx_reader_read(r, a.load_n_threads)) != 0) {
        mm_mapopt_update(&mopt, mi); // this sets the maximum minimizer occurrence; TODO: set a better default in mm_mapopt_init()!
        fprintf(stderr, "[M::%s::%.3f*%.2f] loaded/built the index for %d target sequence(s)\n",
                __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), mi->n_seq);
        lq_collect_occ(mi, a.args[1] == NULL ? 0 : 1, a.args[1] == NULL ? NULL : &a.args[1], &mopt,
                       a.map_n_threads, &m_cnts, &offset_q);
        offset_q   = 0; // all of queries were visited by every chunnk of db.
        mm_idx_destroy(mi);
    }
    
    if(a.args[1] == NULL){
        mm_idx_reader_close(r); // close the index reader
        fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - mm_realtime0, cputime());
        return 0;
    }
    
    // open query file for reading; you may use your favorite FASTA/Q parser
    f = gzopen(a.args[1], "r");
    assert(f);
    ks = kseq_init(f);
    i = 0;
    
    // vars to get minimizers. used for paper revision
    mm128_v minis;
    int h_ret;
    void *km = km_init();
    khiter_t itr;
    uint64_t *b; int32_t n_b;;
    khash_t(64) *h = kh_init(64);
    kh_resize(64, h, 0);
    
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
                uint64_t t;
                t = kh_val(h, itr);
                assert( t == (uint64_t) mv.a[j]);
                //kh_val(h, itr) = kh_val(h, itr) + (uint64_t) mv.a[j];
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
    

    for(j = 0; j < n_b; ++j) printf("minimizer %d cnt: %"PRIu64"\n", j, b[j]);
        
    kfree(km, b);
    km_destroy(km);
    kseq_destroy(ks);
    mm_idx_reader_close(r); // close the index reader
    gzclose(f);
    
    free(m_cnts.a);

    fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - mm_realtime0, cputime());
    return 0;
}

/* Organize functions later neatly*/
typedef struct {
    uint32_t n;
    uint32_t qpos;
    uint32_t seg_id;
    const uint64_t *cr;
} mm_match_t;

struct mm_tbuf_s {
    void *km;
};

static void collect_minimizers(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int n_segs, const int *qlens, const char **seqs, mm128_v *mv)
{
    int i, j, n, sum = 0;
    mv->n = 0;
    for (i = n = 0; i < n_segs; ++i) {
        mm_sketch(km, seqs[i], qlens[i], mi->w, mi->k, i, mi->flag&MM_I_HPC, mv);
        for (j = n; j < mv->n; ++j)
            mv->a[j].y += sum << 1;
        //if (opt->sdust_thres > 0) // mask low-complexity minimizers
        //    mv->n = n + mm_dust_minier(km, mv->n - n, mv->a + n, qlens[i], seqs[i], opt->sdust_thres);
        sum += qlens[i], n = mv->n;
    }
}

static void collect_occurences(const mm_idx_t *mi, mm128_v *mv, m_array *m_cnts)
{
    int i;
    assert(mv->n == m_cnts->n);
    for (i = 0; i < mv->n; ++i) {
        int t;
        mm128_t *p = &mv->a[i];
        mm_idx_get(mi, p->x>>8, &t);
        if(m_cnts->a[i] + t < UINT16_MAX)
            m_cnts->a[i] += t;
    }
}

int lq_collect_occ_frag(const mm_idx_t *mi, int n_segs, const int *qlens, const char **seqs, int *n_regs, mm_tbuf_t *b, const mm_mapopt_t *opt, const char *qname, lq_minimizer_cnt_v *m_cnts, int off_q)
{
    int i, qlen_sum;
    uint32_t hash;
    m_array mc = m_cnts->a[off_q];
    mm128_v mv = {0,0,0};
    km_stat_t kmst;
    
    for (i = 0, qlen_sum = 0; i < n_segs; ++i)
        qlen_sum += qlens[i], n_regs[i] = 0;
    //qlen_sum += qlens[i], n_regs[i] = 0, regs[i] = 0;
    
    if (qlen_sum == 0 || n_segs <= 0 || n_segs > MM_MAX_SEG) return 0;
    
    hash  = qname? __ac_X31_hash_string(qname) : 0;
    hash ^= __ac_Wang_hash(qlen_sum) + __ac_Wang_hash(opt->seed);
    hash  = __ac_Wang_hash(hash);
    
    collect_minimizers(b->km, opt, mi, n_segs, qlens, seqs, &mv);
    // put occurences to m_cnts as db is split by the size specified -I option.
    collect_occurences(mi, &mv, &mc);
    //fprintf(stderr, "debug (count check): seqnum=%d\n", n_segs);

    kfree(b->km, mv.a);
    
    if (b->km) {
        km_stat(b->km, &kmst);
        if (mm_dbg_flag & MM_DBG_PRINT_QNAME)
            fprintf(stderr, "QM\t%s\t%d\tcap=%ld,nCore=%ld,largest=%ld\n", qname, qlen_sum, kmst.capacity, kmst.n_cores, kmst.largest);
        assert(kmst.n_blocks == kmst.n_cores); // otherwise, there is a memory leak
        if (kmst.largest > 1U<<28) {
            km_destroy(b->km);
            b->km = km_init();
        }
    }
    
    return 1;
}

typedef struct {
    int mini_batch_size, n_processed, n_threads, n_fp;
    const mm_mapopt_t *opt;
    mm_bseq_file_t **fp;
    const mm_idx_t *mi;
    kstring_t str;
    int *off_q; // YF memo:
    lq_minimizer_cnt_v *m_cnts; // YF memo:
} pipeline_t;

typedef struct {
    const pipeline_t *p;
    int n_seq, n_frag;
    int off_th_q; // YF memo:
    mm_bseq1_t *seq;
    int *n_reg, *seg_off, *n_seg;
    mm_tbuf_t **buf;
} step_t; //storages for each thred

// worker fucntion for each thread
static void worker_for(void *_data, long i, int tid) // kt_for() callback
{
    step_t *s = (step_t*)_data;
    int qlens[MM_MAX_SEG], j, off = s->seg_off[i], pe_ori = s->p->opt->pe_ori, is_sr = !!(s->p->opt->flag & MM_F_SR);
    //int off_q = off + s->off_th_q; //YF memo: off is the index within a minibatch and s->off_th_q is the offset for the chunk.
    const char *qseqs[MM_MAX_SEG], *quals[MM_MAX_SEG];
    mm_tbuf_t *b  = s->buf[tid];
    //fprintf(stderr, "DEBUG: setting seq2tid. seq=%d, tid=%d\n", off, tid);
    assert(s->n_seg[i] <= MM_MAX_SEG);
    memset(quals, 0, sizeof(char*) * MM_MAX_SEG);
    if (mm_dbg_flag & MM_DBG_PRINT_QNAME)
        fprintf(stderr, "QR\t%s\t%d\t%d\n", s->seq[off].name, tid, s->seq[off].l_seq);
    for (j = 0; j < s->n_seg[i]; ++j) {
        if (s->n_seg[i] == 2 && ((j == 0 && (pe_ori>>1&1)) || (j == 1 && (pe_ori&1))))
            mm_revcomp_bseq(&s->seq[off + j]);
        qlens[j] = s->seq[off + j].l_seq;
        qseqs[j] = s->seq[off + j].seq;
        quals[j] = is_sr? s->seq[off + j].qual : 0;
    }
    if (s->p->opt->flag & MM_F_INDEPEND_SEG) {
        for (j = 0; j < s->n_seg[i]; ++j)
            lq_collect_occ_frag(s->p->mi, 1, &qlens[j], &qseqs[j], &s->n_reg[off+j], b, s->p->opt, s->seq[off+j].name,
                                s->p->m_cnts, off);
    } else {
        lq_collect_occ_frag(s->p->mi, s->n_seg[i], qlens, qseqs, &s->n_reg[off], b, s->p->opt, s->seq[off].name,
                            s->p->m_cnts, off);
        //mm_map_frag(s->p->mi, s->n_seg[i], qlens, qseqs, quals, &s->n_reg[off], &s->reg[off], b, s->p->opt, s->seq[off].name);
    }
}

static void *worker_pipeline(void *shared, int step, void *in)
{
    int i, j, k;
    pipeline_t *p = (pipeline_t*)shared;
    if (step == 0) { // step 0: read sequences
        int with_qual = (!!(p->opt->flag & MM_F_OUT_SAM) && !(p->opt->flag & MM_F_NO_QUAL));
        int frag_mode = (p->n_fp > 1 || !!(p->opt->flag & MM_F_FRAG_MODE));
        step_t *s;
        s = (step_t*)calloc(1, sizeof(step_t));
        if (p->n_fp > 1) s->seq = mm_bseq_read_frag(p->n_fp, p->fp, p->mini_batch_size, with_qual, &s->n_seq);
        else s->seq = mm_bseq_read2(p->fp[0], p->mini_batch_size, with_qual, frag_mode, &s->n_seq);
        if (s->seq) {
            s->p = p;
            for (i = 0; i < s->n_seq; ++i)
                s->seq[i].rid = p->n_processed++;
            s->buf  = (mm_tbuf_t**)calloc(p->n_threads, sizeof(mm_tbuf_t*));
            for (i = 0; i < p->n_threads; ++i){
                s->buf[i] = mm_tbuf_init();
            }
            s->n_reg = (int*)calloc(3 * s->n_seq, sizeof(int));
            s->seg_off = s->n_reg + s->n_seq; // seg_off and n_seg are allocated together with n_reg
            s->n_seg = s->seg_off + s->n_seq;
            for (i = 1, j = 0; i <= s->n_seq; ++i)
                if (i == s->n_seq || !frag_mode || !mm_qname_same(s->seq[i-1].name, s->seq[i].name)) {
                    s->n_seg[s->n_frag] = i - j;
                    s->seg_off[s->n_frag++] = j;
                    j = i;
                }
            //printf("DEBUG: off_q = %d, n_seq=%d, step offset=%d\n", *p->off_q, s->n_seq, s->off_th_q);
            return s;
        } else {
            //fprintf(stderr, "DEBUG: no more seq off_q = %d\n", *p->off_q);
            free(s);
        }
    } else if (step == 1) { // step 1: map
        kt_for(p->n_threads, worker_for, in, ((step_t*)in)->n_frag);
        return in;
    } else if (step == 2) { // step 2: output
        step_t *s = (step_t*)in;
        for (i = 0; i < p->n_threads; ++i){
            mm_tbuf_destroy(s->buf[i]);
        }
        free(s->buf);
        for (k = 0; k < s->n_frag; ++k) {
            int seg_st = s->seg_off[k], seg_en = s->seg_off[k] + s->n_seg[k];
            for (i = seg_st; i < seg_en; ++i) {
                //for (j = 0; j < s->n_reg[i]; ++j) free(s->reg[i][j].p);
                //free(s->reg[i]);
                free(s->seq[i].seq); free(s->seq[i].name);
                if (s->seq[i].qual) free(s->seq[i].qual);
            }
        }
        free(s->n_reg); free(s->seq); // seg_off and n_seg were allocated with reg; no memory leak here
        fprintf(stderr, "[M::%s::%.3f*%.2f] processed %d sequences. (Peak RSS: %.3f GB) \n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), s->n_seq, peakrss() / 1024.0 / 1024.0 / 1024.0);
        free(s);
    }
    return 0;
}

int lq_collect_occ(const mm_idx_t *idx, int n_segs, const char **fn, const mm_mapopt_t *opt,
                   int n_threads, lq_minimizer_cnt_v *m_cnts, int *off)
{
    int i, j, pl_threads;
    pipeline_t pl;
    if (n_segs < 1) return -1;
    memset(&pl, 0, sizeof(pipeline_t));
    pl.n_fp = n_segs;
    pl.fp = (mm_bseq_file_t**)calloc(n_segs, sizeof(mm_bseq_file_t*));
    for (i = 0; i < n_segs; ++i) {
        pl.fp[i] = mm_bseq_open(fn[i]);
        if (pl.fp[i] == 0) {
            if (mm_verbose >= 1)
                fprintf(stderr, "ERROR: failed to open file '%s'\n", fn[i]);
            for (j = 0; j < i; ++j)
                mm_bseq_close(pl.fp[j]);
            free(pl.fp);
            return -1;
        }
    }
    pl.opt = opt, pl.mi = idx;
    pl.n_threads = n_threads > 1? n_threads : 1;
    pl.mini_batch_size = opt->mini_batch_size;
    pl_threads = n_threads == 1? 1 : (opt->flag&MM_F_2_IO_THREADS)? 3 : 2;
    pl.off_q   = off; //YF memo
    pl.m_cnts = m_cnts; // YF:
    kt_pipeline(pl_threads, worker_pipeline, &pl, 3);
    free(pl.str.s);
    for (i = 0; i < n_segs; ++i)
        mm_bseq_close(pl.fp[i]);
    free(pl.fp);
    return 0;
}

