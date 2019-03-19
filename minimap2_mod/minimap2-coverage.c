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
#include "lqfuncs.h"
#include "lqutils.h"
#include "kthread.h"
#include "kvec.h"
#include "minimap2-coverage.h"
#include "ksort.h"
#include <unistd.h>

#define sort_key_32(x) (x)
KRADIX_SORT_INIT(32, uint32_t, sort_key_32, 4)

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

static void chain_post(const mm_mapopt_t *opt, int max_chain_gap_ref, const mm_idx_t *mi, void *km, int qlen, int n_segs, const int *qlens, int *n_regs, mm_reg1_t *regs, mm128_t *a)
{
    if (!(opt->flag & (MM_F_AVA|LQ_F_AVA))) { // don't choose primary mapping(s) for read overlap
        mm_set_parent(km, opt->mask_level, *n_regs, regs, opt->a * 2 + opt->b);
        if (n_segs <= 1) mm_select_sub(km, opt->pri_ratio, mi->k*2, opt->best_n, n_regs, regs);
        else mm_select_sub_multi(km, opt->pri_ratio, 0.2f, 0.7f, max_chain_gap_ref, mi->k*2, opt->best_n, n_segs, qlens, n_regs, regs);
        if (!(opt->flag & MM_F_SPLICE) && !(opt->flag & MM_F_SR) && !(opt->flag & MM_F_NO_LJOIN))
            mm_join_long(km, opt, qlen, n_regs, regs, a);
    }
}

/*
static mm_reg1_t *align_regs(const mm_mapopt_t *opt, const mm_idx_t *mi, void *km, int qlen, const char *seq, const char *qual, int *n_regs, mm_reg1_t *regs, mm128_t *a)
{
    if (!(opt->flag & MM_F_CIGAR)) return regs;
    regs = mm_align_skeleton(km, opt, mi, qlen, seq, n_regs, regs, a); // this calls mm_filter_regs()
    if (!(opt->flag & (MM_F_AVA|LQ_F_AVA))) {
        mm_set_parent(km, opt->mask_level, *n_regs, regs, opt->a * 2 + opt->b);
        mm_select_sub(km, opt->pri_ratio, mi->k*2, opt->best_n, n_regs, regs);
        mm_set_sam_pri(*n_regs, regs);
    }
    return regs;
}
*/

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

static mm128_t *collect_seed_hits(void *km, const mm_mapopt_t *opt, int max_occ, const mm_idx_t *mi, const char *qname, const mm128_v *mv, int qlen, int64_t *n_a, int *rep_len,
                                  int *n_mini_pos, uint64_t **mini_pos)
{
    int rep_st = 0, rep_en = 0, i;
    mm_match_t *m;
    mm128_t *a;
    
    *n_mini_pos = 0;
    *mini_pos = (uint64_t*)kmalloc(km, mv->n * sizeof(uint64_t));
    m = (mm_match_t*)kmalloc(km, mv->n * sizeof(mm_match_t));
    for (i = 0; i < mv->n; ++i) {
        int t;
        mm128_t *p = &mv->a[i];
        m[i].qpos = (uint32_t)p->y;
        m[i].cr = mm_idx_get(mi, p->x>>8, &t);
        m[i].n = t; // YF memo: occurenace of i-th minimizer in the db?
        m[i].seg_id = p->y >> 32;
    }
    for (i = 0, *n_a = 0; i < mv->n; ++i) // find the length of a[]
        if (m[i].n < max_occ) *n_a += m[i].n;
    a = (mm128_t*)kmalloc(km, *n_a * sizeof(mm128_t));
    for (i = *rep_len = 0, *n_a = 0; i < mv->n; ++i) {
        mm128_t *p = &mv->a[i];
        mm_match_t *q = &m[i];
        const uint64_t *r = q->cr;
        int k, q_span = p->x & 0xff, is_tandem = 0;
        if (q->n >= max_occ) {
            int en = (q->qpos>>1) + 1, st = en - q_span;
            if (st > rep_en) {
                *rep_len += rep_en - rep_st;
                rep_st = st, rep_en = en;
            } else rep_en = en;
            continue;
        }
        (*mini_pos)[(*n_mini_pos)++] = (uint64_t)q_span<<32 | q->qpos>>1;
        if (i > 0 && p->x>>8 == mv->a[i - 1].x>>8) is_tandem = 1;
        if (i < mv->n - 1 && p->x>>8 == mv->a[i + 1].x>>8) is_tandem = 1;
        for (k = 0; k < q->n; ++k) {
            int32_t rpos = (uint32_t)r[k] >> 1;
            mm128_t *p;
            if (qname && (opt->flag&(MM_F_NO_SELF|MM_F_AVA))) {
                const char *tname = mi->seq[r[k]>>32].name;
                //fprintf(stderr, "DEBUG tname = %s, rid = %llu\n", tname, r[k] >> 32);
                int cmp;
                cmp = strcmp(qname, tname);
                if ((opt->flag&MM_F_NO_SELF) && cmp == 0 && rpos == (q->qpos>>1)) // avoid the diagonal
                    continue;
                if ((opt->flag&MM_F_AVA) && cmp > 0) // all-vs-all mode: map once
                    continue;
            }
            p = &a[(*n_a)++];
            if ((r[k]&1) == (q->qpos&1)) { // forward strand
                p->x = (r[k]&0xffffffff00000000ULL) | rpos;
                p->y = (uint64_t)q_span << 32 | q->qpos >> 1;
            } else { // reverse strand
                p->x = 1ULL<<63 | (r[k]&0xffffffff00000000ULL) | rpos;
                p->y = (uint64_t)q_span << 32 | (qlen - ((q->qpos>>1) + 1 - q_span) - 1);
            }
            p->y |= (uint64_t)q->seg_id << MM_SEED_SEG_SHIFT;
            if (is_tandem) p->y |= MM_SEED_TANDEM;
        }
    }
    *rep_len += rep_en - rep_st;
    kfree(km, m);
    return a;
}

int lq_map_frag(const mm_idx_t *mi, int n_segs, const int *qlens, const char **seqs, const char **quals, int *n_regs, mm_reg1_t **regs, mm_tbuf_t *b, const mm_mapopt_t *opt, const lq_fltopt_t *fopt, int min_score_med_good, const char *qname, uint64_t *lambdas, uint64_t *lambdas2, lq_subcoords_v **ovlp_cors, lq_minimizer_cnt_v **m_cnts, float *ks, uint32_t off_db, uint32_t off_q)
{
    int i, rep_len, qlen_sum, n_regs0, n_mini_pos;
    
    int max_chain_gap_qry, max_chain_gap_ref, is_splice = !!(opt->flag & MM_F_SPLICE), is_sr = !!(opt->flag & MM_F_SR);
    uint32_t hash;
    int64_t n_a;
    uint64_t *u, *mini_pos;

    // YF memo:
    if(ovlp_cors[off_q] == NULL){
        ovlp_cors[off_q] = (lq_subcoords_v *)malloc(sizeof(lq_subcoords_v));
        kv_init(*ovlp_cors[off_q]);
    }
    lq_subcoords_v *vc = ovlp_cors[off_q];
    
    // temporary holder
    uint64_t lambda_temp;
    uint64_t lambda2_temp;
    kvec_t(lq_tchain) v = {0,0,0};
    
    mm128_t *a;
    mm128_v mv = {0,0,0};
    mm_reg1_t *regs0;
    km_stat_t kmst;
    
    for (i = 0, qlen_sum = 0; i < n_segs; ++i)
        qlen_sum += qlens[i], n_regs[i] = 0;
    //qlen_sum += qlens[i], n_regs[i] = 0, regs[i] = 0;
    
    if (qlen_sum == 0 || n_segs <= 0 || n_segs > MM_MAX_SEG) return 0;
    
    hash  = qname? __ac_X31_hash_string(qname) : 0;
    hash ^= __ac_Wang_hash(qlen_sum) + __ac_Wang_hash(opt->seed);
    hash  = __ac_Wang_hash(hash);
    
    collect_minimizers(b->km, opt, mi, n_segs, qlens, seqs, &mv);
    a = collect_seed_hits(b->km, opt, opt->mid_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);
    radix_sort_128x(a, a + n_a);
    
    if(m_cnts[off_q] == NULL){
        m_cnts[off_q] = (lq_minimizer_cnt_v *)malloc(sizeof(lq_minimizer_cnt_v));
        kv_init(*m_cnts[off_q]);
        // YF: counter reset for the first batch
        for(i=0; i<mv.n; i++)
            kv_push(uint32_t, 0, *m_cnts[off_q], 0);
    }
    lq_minimizer_cnt_v *mc = m_cnts[off_q];
    
    // set max chaining gap on the query and the reference sequence
    if (is_sr)
        max_chain_gap_qry = qlen_sum > opt->max_gap? qlen_sum : opt->max_gap;
    else max_chain_gap_qry = opt->max_gap;
    if (opt->max_gap_ref > 0) {
        max_chain_gap_ref = opt->max_gap_ref; // always honor mm_mapopt_t::max_gap_ref if set
    } else if (opt->max_frag_len > 0) {
        max_chain_gap_ref = opt->max_frag_len - qlen_sum;
        if (max_chain_gap_ref < opt->max_gap) max_chain_gap_ref = opt->max_gap;
    } else max_chain_gap_ref = opt->max_gap;
    
    lambda_temp  = lambdas[off_q];
    lambda2_temp = lambdas2[off_q];
    a = lq_chain_dp(mi, max_chain_gap_ref, max_chain_gap_qry, fopt, opt->bw, opt->max_chain_skip, opt->min_cnt, opt->min_chain_score, min_score_med_good, is_splice, n_segs, n_a, a, &n_regs0, &u, b->km, qlen_sum, &lambdas[off_q], &lambdas2[off_q], &v, vc);
    
    if (opt->max_occ > opt->mid_occ && rep_len > 0) {
        int rechain = 0;
        if (n_regs0 > 0) { // test if the best chain has all the segments
            int n_chained_segs = 1, max = 0, max_i = -1, max_off = -1, off = 0;
            for (i = 0; i < n_regs0; ++i) { // find the best chain
                if (max < u[i]>>32) max = u[i]>>32, max_i = i, max_off = off;
                off += (uint32_t)u[i];
            }
            for (i = 1; i < (uint32_t)u[max_i]; ++i) // count the number of segments in the best chain
                if ((a[max_off+i].y&MM_SEED_SEG_MASK) != (a[max_off+i-1].y&MM_SEED_SEG_MASK))
                    ++n_chained_segs;
            if (n_chained_segs < n_segs)
                rechain = 1;
        } else rechain = 1;
        if (rechain) { // redo chaining with a higher max_occ threshold
            kfree(b->km, a);
            kfree(b->km, u);
            kfree(b->km, mini_pos);
            a = collect_seed_hits(b->km, opt, opt->max_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);
            radix_sort_128x(a, a + n_a);
            lambdas[off_q]  = lambda_temp;  // YF memo: recomp
            lambdas2[off_q] = lambda2_temp; // YF memo: recomp
            kv_init(v); // YF memo: recomp
            kv_init(*vc);
            a = lq_chain_dp(mi, max_chain_gap_ref, max_chain_gap_qry, fopt, opt->bw, opt->max_chain_skip, opt->min_cnt, opt->min_chain_score, min_score_med_good, is_splice, n_segs, n_a, a, &n_regs0, &u, b->km, qlen_sum, &lambdas[off_q], &lambdas2[off_q], &v, vc);
        }
    }
    
    if(opt->flag&MM_F_AVA){
        for(i = 0; i < v.n; i++){
            lambdas[v.a[i].rid + off_db] += v.a[i].length;
        }
    }

    // followed relevant code from miniasm
    // the rest part is written in compute_reliable_region()
    radix_sort_32(vc->a, vc->a + vc->n);

    kv_destroy(v);
    
    //lq_cnt_exact_match(min_score_med_good, qlen_sum, n_regs0, u, a, n_mini_pos, mini_pos, mc);
    
    regs0 = mm_gen_regs(b->km, hash, qlen_sum, n_regs0, u, a);
    
    chain_post(opt, max_chain_gap_ref, mi, b->km, qlen_sum, n_segs, qlens, &n_regs0, regs0, a);
    if (!is_sr) lq_cnt_match(mi, qlen_sum, lambdas[off_q], n_regs0, regs0, a, n_mini_pos, mini_pos, mc, min_score_med_good, &ks[off_q], fopt);

    /*
    if (n_segs == 1) { // uni-segment
        regs0 = align_regs(opt, mi, b->km, qlens[0], seqs[0], quals? quals[0] : 0, &n_regs0, regs0, a);
        mm_set_mapq(n_regs0, regs0, opt->min_chain_score, opt->a, rep_len, is_sr);
        n_regs[0] = n_regs0, regs[0] = regs0;
    } else { // multi-segment
        //mm_seg_t *seg;
        //seg = mm_seg_gen(b->km, hash, n_segs, qlens, n_regs0, regs0, n_regs, regs, a); // split fragment chain to separate segment chains
        //free(regs0);
        //for (i = 0; i < n_segs; ++i) {
        //    mm_set_parent(b->km, opt->mask_level, n_regs[i], regs[i], opt->a * 2 + opt->b); // update mm_reg1_t::parent
        //    regs[i] = align_regs(opt, mi, b->km, qlens[i], seqs[i], quals? quals[i] : 0, &n_regs[i], regs[i], seg[i].a);
        //    mm_set_mapq(n_regs[i], regs[i], opt->min_chain_score, opt->a, rep_len, is_sr);
        //}
        //mm_seg_free(b->km, n_segs, seg);
        //if (n_segs == 2 && opt->pe_ori >= 0 && (opt->flag&MM_F_CIGAR))
        //   mm_pair(b->km, max_chain_gap_ref, opt->pe_bonus, opt->a * 2 + opt->b, opt->a, qlens, n_regs, regs); // pairing
    }
    */
    
    kfree(b->km, mv.a);
    kfree(b->km, a);
    kfree(b->km, u);
    kfree(b->km, mini_pos);
    
    if (b->km) {
        km_stat(b->km, &kmst);
        if (mm_dbg_flag & MM_DBG_PRINT_QNAME)
            fprintf(stderr, "QM\t%s\t%d\tcap=%ld,nCore=%ld,largest=%ld\n", qname, qlen_sum, kmst.capacity, kmst.n_cores, kmst.largest);
        assert(kmst.n_blocks == kmst.n_cores); // otherwise, there is a memory leak
        if (kmst.largest > 1U<<28) {
            km_destroy(b->km);
            b->km = km_init();
            fprintf(stderr, "DEBUG: km was reinitialized at off_db: %d, off_q: %d\n", off_db, off_q);
        }
    }
    
    return 1;
}

//**************************
//* Multi-threaded mapping *
//**************************

typedef struct {
    int mini_batch_size, n_processed, n_threads, n_fp;
    const mm_mapopt_t *opt;
    mm_bseq_file_t **fp;
    const mm_idx_t *mi;
    kstring_t str;
    const lq_fltopt_t *fopt; // YF memo:
    uint32_t off_db; //YF memo: database offset. This is updated at most outer loop.
    uint32_t *off_q; // YF memo:
    uint32_t min_score_med_good; //YF memo:
    uint64_t *lambdas; //YF memo:
    uint64_t *lambdas2; //YF memo:
    float *avg_ks; //YF memo:
    lq_subcoords_v **ovlp_cors; //YF memo:
    lq_minimizer_cnt_v **m_cnts; // YF memo:
} pipeline_t;

typedef struct {
    const pipeline_t *p;
    int n_seq, n_frag;
    uint32_t offset_q; // YF memo:
    mm_bseq1_t *seq;
    int *n_reg, *seg_off, *n_seg;
    mm_reg1_t **reg;
    mm_tbuf_t **buf;
} step_t;

static void worker_for(void *_data, long i, int tid) // kt_for() callback
{
    step_t *s = (step_t*)_data;
    int qlens[MM_MAX_SEG], j, off = s->seg_off[i], pe_ori = s->p->opt->pe_ori, is_sr = !!(s->p->opt->flag & MM_F_SR);
    uint32_t off_q = off + s->offset_q; // off is the index within the chunk and s->offset_q is the offset for the chunk.
    const char *qseqs[MM_MAX_SEG], *quals[MM_MAX_SEG];
    mm_tbuf_t *b = s->buf[tid];
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
            lq_map_frag(s->p->mi, 1, &qlens[j], &qseqs[j], &quals[j], &s->n_reg[off+j], &s->reg[off+j], b, s->p->opt, s->p->fopt, s->p->min_score_med_good, s->seq[off+j].name, s->p->lambdas, s->p->lambdas2,s->p->ovlp_cors, s->p->m_cnts, s->p->avg_ks, s->p->off_db, off_q);
    } else {
        lq_map_frag(s->p->mi, s->n_seg[i], qlens, qseqs, quals, &s->n_reg[off], &s->reg[off], b, s->p->opt, s->p->fopt, s->p->min_score_med_good, s->seq[off].name, s->p->lambdas, s->p->lambdas2, s->p->ovlp_cors, s->p->m_cnts, s->p->avg_ks, s->p->off_db, off_q);
    }
    for (j = 0; j < s->n_seg[i]; ++j) // flip the query strand and coordinate to the original read strand
        if (s->n_seg[i] == 2 && ((j == 0 && (pe_ori>>1&1)) || (j == 1 && (pe_ori&1)))) {
            int k, t;
            mm_revcomp_bseq(&s->seq[off + j]);
            for (k = 0; k < s->n_reg[off + j]; ++k) {
                mm_reg1_t *r = &s->reg[off + j][k];
                t = r->qs;
                r->qs = qlens[j] - r->qe;
                r->qe = qlens[j] - t;
                r->rev = !r->rev;
            }
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
            s->buf = (mm_tbuf_t**)calloc(p->n_threads, sizeof(mm_tbuf_t*));
            for (i = 0; i < p->n_threads; ++i)
                s->buf[i] = mm_tbuf_init();
            s->n_reg = (int*)calloc(3 * s->n_seq, sizeof(int));
            s->seg_off = s->n_reg + s->n_seq; // seg_off and n_seg are allocated together with n_reg
            s->n_seg = s->seg_off + s->n_seq;
            s->reg = (mm_reg1_t**)calloc(s->n_seq, sizeof(mm_reg1_t*));
            for (i = 1, j = 0; i <= s->n_seq; ++i)
                if (i == s->n_seq || !frag_mode || !mm_qname_same(s->seq[i-1].name, s->seq[i].name)) {
                    s->n_seg[s->n_frag] = i - j;
                    s->seg_off[s->n_frag++] = j;
                    j = i;
                }
            s->offset_q  = *p->off_q;
            *(p->off_q) += s->n_seq; // add offset of query.
            //printf("DEBUG: off_q = %d, n_seq=%d, step offset=%d\n", *p->off_q, s->n_seq, s->offset_q);
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
        for (i = 0; i < p->n_threads; ++i) mm_tbuf_destroy(s->buf[i]);
        free(s->buf);
        for (k = 0; k < s->n_frag; ++k) {
            int seg_st = s->seg_off[k], seg_en = s->seg_off[k] + s->n_seg[k];
	    for (i = seg_st; i < seg_en; ++i) {
                for (j = 0; j < s->n_reg[i]; ++j) free(s->reg[i][j].p);
                free(s->reg[i]);
                free(s->seq[i].seq); free(s->seq[i].name);
                if (s->seq[i].qual) free(s->seq[i].qual);
            }
        }
        free(s->reg); free(s->n_reg); free(s->seq); // seg_off and n_seg were allocated with reg; no memory leak here
        fprintf(stderr, "[M::%s::%.3f*%.2f] mapped %d sequences. (Peak RSS: %.3f GB) \n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), s->n_seq, peakrss() / 1024.0 / 1024.0 / 1024.0);
        free(s);
    }
    return 0;
}

// kovlp_cors is vec_t(lq_subcoords_t)*
int lq_map_file_frag(const mm_idx_t *idx, int n_segs, const char **fn, const mm_mapopt_t *opt, const lq_fltopt_t *fopt, int min_score_med, int min_score_good, int n_threads, uint64_t *lambdas, uint64_t *lambdas2, lq_subcoords_v **ovlp_cors, lq_minimizer_cnt_v **m_cnts, float *ks, uint32_t off_db, uint32_t *off)
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
    pl.off_db  = off_db; // YF memo
    pl.off_q   = off; //YF memo
    pl.lambdas = lambdas; // YF memo:
    pl.lambdas2 = lambdas2; // YF memo: for primary chain
    pl.ovlp_cors = ovlp_cors; // YF memo:
    pl.m_cnts = m_cnts; // YF:
    pl.avg_ks = ks; // YF
    pl.min_score_med_good = (min_score_med << 16) | min_score_good; // YF memo
    pl.fopt = fopt;
    kt_pipeline(pl_threads, worker_pipeline, &pl, 3);
    free(pl.str.s);
    for (i = 0; i < n_segs; ++i)
        mm_bseq_close(pl.fp[i]);
    free(pl.fp);
    return 0;
}

int lq_map_file(const mm_idx_t *idx, const char *fn, const mm_mapopt_t *opt, const lq_fltopt_t *fopt, int min_score_med, int min_score_good, int n_threads, uint64_t *lambdas, uint64_t *lambdas2, lq_subcoords_v **ovlp_cors, lq_minimizer_cnt_v **m_cnts, float *ks, uint32_t off_db, uint32_t *off)
{
    int fnum = (fn == NULL) ? 0 : 1;
    return lq_map_file_frag(idx, fnum, &fn, opt, fopt, min_score_med, min_score_good, n_threads, lambdas, lambdas2, ovlp_cors, m_cnts, ks, off_db, off);
}

// so dirty. need more revision
void compute_reliable_region(lq_subcoords_v *vc, uint32_t min_cov, lq_region_v *coords, lq_region_v *mcoords, lq_region_v *gcoords)
{
    int j;
    uint32_t start, cov;
    uint32_t med_start, med_cov;
    uint32_t good_start, good_cov;
    start = cov = med_start = med_cov = good_start = good_cov = 0;
    
    for (j = 0, cov = 0; j < vc->n; ++j) {
        uint32_t old_cov      = cov;
        uint32_t old_med_cov  = med_cov;
        uint32_t old_good_cov = good_cov;
        if (vc->a[j]&1){
            --cov;
            if(vc->a[j]&2)
                --med_cov;
            if(vc->a[j]&4)
                --good_cov;
        }
        else {
            ++cov;
            if(vc->a[j]&2)
                ++med_cov;
            if(vc->a[j]&4)
                ++good_cov;
        }
        if (old_cov < min_cov && cov >= min_cov) {
            start = vc->a[j]>>3;
            // just in case. cov and med_cov and good_cov is not 100% exclusive.
            if(old_med_cov < min_cov && med_cov >= min_cov){
                med_start = vc->a[j]>>3;
                if (old_good_cov < min_cov && good_cov >= min_cov)
                    good_start = vc->a[j]>>3;
            }
        } else if (old_cov >= min_cov && cov < min_cov) {
            uint32_t len = (vc->a[j]>>3) - start;
            if (len > 0){
                lq_subcoords_t c;
                c.start = start; c.end = vc->a[j]>>3;
                kv_push(lq_subcoords_t, 0, *coords, c);
            }
            // just in case. cov and med_cov and good_cov is not 100% exclusive.
            if (old_med_cov >= min_cov && med_cov < min_cov){
                uint32_t mlen = (vc->a[j]>>3) - med_start;
                if(mlen > 0){
                    lq_subcoords_t mc;
                    mc.start = med_start; mc.end = vc->a[j]>>3;
                    kv_push(lq_subcoords_t, 0, *mcoords, mc);
                }
                if(old_good_cov >= min_cov && good_cov < min_cov){
                    uint32_t glen = (vc->a[j]>>3) - good_start;
                    if(glen > 0){
                        lq_subcoords_t gc;
                        gc.start = good_start; gc.end = vc->a[j]>>3;
                        kv_push(lq_subcoords_t, 0, *gcoords, gc);
                    }
                }
            }
        } else if (old_med_cov < min_cov && med_cov >= min_cov) {
            med_start = vc->a[j]>>3;
            // just in case. med_cov and good_cov is not 100% exclusive.
            if (old_good_cov < min_cov && good_cov >= min_cov)
                good_start = vc->a[j]>>3;
        } else if(old_med_cov >= min_cov && med_cov < min_cov) {
            uint32_t mlen = (vc->a[j]>>3) - med_start;
            if(mlen > 0){
                lq_subcoords_t mc;
                mc.start = med_start; mc.end = vc->a[j]>>3;
                kv_push(lq_subcoords_t, 0, *mcoords, mc);
            }
            // just in case. med_cov and good_cov is not 100% exclusive.
            if(old_good_cov >= min_cov && good_cov < min_cov){
                uint32_t glen = (vc->a[j]>>3) - good_start;
                if(glen > 0){
                    lq_subcoords_t gc;
                    gc.start = good_start; gc.end = vc->a[j]>>3;
                    kv_push(lq_subcoords_t, 0, *gcoords, gc);
                }
            }
        } else if (old_good_cov < min_cov && good_cov >= min_cov) {
            good_start = vc->a[j]>>3;
        } else if (old_good_cov >= min_cov && good_cov < min_cov){
            uint32_t glen = (vc->a[j]>>3) - good_start;
            if(glen > 0){
                lq_subcoords_t gc;
                gc.start = good_start; gc.end = vc->a[j]>>3;
                kv_push(lq_subcoords_t, 0, *gcoords, gc);
            }
        }
    }
}


/* End of messy functions*/

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

const char *argp_program_version = "minimap2-coverage 0.2 (forked from 2.6-r639-dirty)";
const char *argp_program_bug_address = "<yoshinori.fukasawa@kaust.edu.sa>";
//static char usage[] = "minimap2-coverage [options] <target.seqs> <query.seqs>";

struct arguments
{
    short h_flag, ava_flag, avs_flag, filter_flag;
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
    int i, j, k, l;
    uint32_t offset_db = 0; // offset for db index
    uint32_t offset_q  = 0;

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
        a.min_coverage = 3;
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
    fprintf(stderr, "Homo-polymer compression: %d, Filtering: %d\n", a.h_flag, a.filter_flag);
    fprintf(stderr, "max-overhang %d, min-overlaplen %d, min-overapratio %.2f\n", fopt.max_overhang, fopt.min_ovlp, fopt.min_ratio);
    fprintf(stderr, "num of threads %d, num of query seqs %d\n===\n", a.map_n_threads, a.num_subsetseq);
    
    // coverage counter
    // YF memo: each step is protected by mutex (see ktp_worker in the kthread.c)
    uint64_t *lambdas; lambdas = (uint64_t*) calloc(a.num_subsetseq, sizeof(uint64_t));
    uint64_t *lambdas2; lambdas2 = (uint64_t*) calloc(a.num_subsetseq, sizeof(uint64_t)); // primary chain coverage
    lq_subcoords_v **ovlp_coords; ovlp_coords = (lq_subcoords_v**) malloc(a.num_subsetseq * sizeof(lq_subcoords_v*));
    lq_minimizer_cnt_v **m_cnts; m_cnts = (lq_minimizer_cnt_v**) malloc(a.num_subsetseq * sizeof(lq_minimizer_cnt_v*));
    float *avg_ks; avg_ks = (float*) calloc(a.num_subsetseq, sizeof(float));
    
    // index can be read several times due to the memory size specified -I option. see iopt.batch_size
    // That's why this has to be in a while loop.
    
    while ((mi = mm_idx_reader_read(r, a.load_n_threads)) != 0) {
        mm_mapopt_update(&mopt, mi); // this sets the maximum minimizer occurrence; TODO: set a better default in mm_mapopt_init()!
        fprintf(stderr, "[M::%s::%.3f*%.2f] loaded/built the index for %d target sequence(s)\n",
                __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), mi->n_seq);
        lq_map_file(mi, a.args[1] == NULL ? NULL : a.args[1], &mopt, &fopt, a.min_score_med, a.min_score_good, a.map_n_threads, lambdas, lambdas2, ovlp_coords, m_cnts, avg_ks, offset_db, &offset_q);
        offset_db += mi->n_seq;
        offset_q   = 0; // all of queries were visited by every chunnk of db.
        mm_idx_destroy(mi);
        
        //fprintf(stderr, "[M::%s::%.3f*%.2f] Releasing freed memory...\n",
        //        __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0));
        //MallocExtension_ReleaseFreeMemory(); // temporal solution.
        //fprintf(stderr, "[M::%s::%.3f*%.2f] Releasing freed memory is done.\n",
        //        __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0));
    }
    
    if(a.args[1] == NULL){
        free(lambdas);
        free(lambdas2);
        free(ovlp_coords);
        free(m_cnts);
        free(avg_ks);
        mm_idx_reader_close(r); // close the index reader
        fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - mm_realtime0, cputime());
        return 0;
    }
    
    // open query file for reading; you may use your favorite FASTA/Q parser
    gzFile f = gzopen(a.args[1], "r");
    assert(f);
    kseq_t *ks = kseq_init(f);
    i = 0;
    while (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence
        lq_region_v regs  = {0,0,0};
        lq_region_v mregs = {0,0,0};
        lq_region_v gregs = {0,0,0};
        lq_subcoords_v* v = ovlp_coords[i];
        
        lq_minimizer_cnt_v* mv = m_cnts[i];
        int32_t n_match = 0; double div = 0.0;
        uint32_t sum = 0.0;
        //fprintf(stderr, "debug-match-cnt %s: ", ks->name.s);
        for(j=0; j<mv->n; j++){
            sum += mv->a[j];
        }
        sum /= mv->n;
        for(j=0; j<mv->n; j++){
            if(mv->a[j] > sum) n_match++;
            //fprintf(stderr, "%d, ", mv->a[j]);
        }
        //fprintf(stderr, "ave. %d\n", sum/mv->n);
        div = n_match > 0 ? logf((float)mv->n / n_match) / avg_ks[i] : 1.0;
        
        compute_reliable_region(v, a.min_coverage, &regs, &mregs, &gregs);
        
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
                if (gregs.n > 0){
                    char *gcoords; gcoords = (char *) malloc(sizeof(char) * (66 *v->n)); gcoords[0] = '\0';
                    for(l=0; l<gregs.n; l++){
                        if(l>0)
                            strcat(gcoords, ",");
                        char gc[66] = {'\0'}; sprintf(gc, "%d-%d", gregs.a[l].start, gregs.a[l].end);
                        strcat(gcoords, gc);
                    }
                    if(a.filter_flag){
                        printf("%s\t%d\t%"PRIu64"\t%s\t%s\t%s\t%.3f\t%.3f\t%.3f\t0.0\n", ks->name.s, ks->seq.l, lambdas[i],
                               coords, mcoords, gcoords, (double) tot/ks->seq.l, meanQ(ks->qual.s, ks->qual.l), div);
                    } else {
                        printf("%s\t%d\t%"PRIu64"\t%s\t%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\n", ks->name.s, ks->seq.l, lambdas[i],
                               coords, mcoords, gcoords, (double) lambdas[i]/tot, meanQ(ks->qual.s, ks->qual.l), div, (double) lambdas2[i]/tot);
                    }
                } else {
                    if(a.filter_flag){
                        printf("%s\t%d\t%"PRIu64"\t%s\t%s\t0\t%.3f\t%.3f\t%.3f\t0.0\n", ks->name.s, ks->seq.l, lambdas[i],
                               coords, mcoords, (double) tot/ks->seq.l, meanQ(ks->qual.s, ks->qual.l), div);
                    } else {
                        printf("%s\t%d\t%"PRIu64"\t%s\t%s\t0\t%.3f\t%.3f\t%.3f\t%.3f\n", ks->name.s, ks->seq.l, lambdas[i],
                               coords, mcoords, (double) lambdas[i]/tot, meanQ(ks->qual.s, ks->qual.l), div, (double) lambdas2[i]/tot);
                    }
                }
            } else {
                if(a.filter_flag){
                    printf("%s\t%d\t%"PRIu64"\t%s\t0\t0\t%.3f\t%.3f\t%.3f\t0.0\n", ks->name.s, ks->seq.l, lambdas[i],
                           coords, (double) tot/ks->seq.l, meanQ(ks->qual.s, ks->qual.l), div);
                } else {
                    printf("%s\t%d\t%"PRIu64"\t%s\t0\t0\t%.3f\t%.3f\t%.3f\t%.3f\n", ks->name.s, ks->seq.l, lambdas[i],
                           coords, (double) lambdas[i]/tot, meanQ(ks->qual.s, ks->qual.l), div, (double) lambdas2[i]/tot);
                }
            }
        } else {
            printf("%s\t%d\t%"PRIu64"\t0\t0\t0\t0.0\t%.3f\t%.3f\t0.0\n", ks->name.s, ks->seq.l, lambdas[i], meanQ(ks->qual.s, ks->qual.l), div);
        }
        if(regs.n > 0)
            kv_destroy(regs);
        if(mregs.n > 0)
            kv_destroy(mregs);
        if(gregs.n > 0)
            kv_destroy(gregs);
        if(v->n > 0)
            kv_destroy(*v);
        free(v);
        
        if(mv->n>0)
            kv_destroy(*mv);
        free(mv);
        
        i++;
    }
    
    free(lambdas);
    free(lambdas2);
    free(ovlp_coords);
    free(m_cnts);
    free(avg_ks);

    kseq_destroy(ks);
        
    mm_idx_reader_close(r); // close the index reader
    gzclose(f);
    
    fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - mm_realtime0, cputime());
    return 0;
}
