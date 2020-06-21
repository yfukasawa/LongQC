#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "kthread.h"
#include "kvec.h"
#include "kalloc.h"
#include "sdust.h"
#include "mmpriv.h"
#include "bseq.h"
#include "khash.h"
#include "minimap2-coverage.h"

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

void filter_redundant_coords(lq_subcoords_v *v, lq_subcoords_v *cv, void *km, uint32_t min_cov)
{
    int i, j;
    lq_subcoords_v mcoords = {0,0,0};
    uint32_t med_start, med_cov;
    med_start = med_cov = 0;
    
    if(!cv || cv->n == 0) return;
    
    kvec_t(uint32_t) vc = {0,0,0};
    for(i=0; i<cv->n;++i){
        kv_push(uint32_t, km, vc, cv->a[i].start);
        kv_push(uint32_t, km, vc, cv->a[i].end);
    }
    radix_sort_32(vc.a, vc.a + vc.n);
   
    // YF memo: reads having a higher score than GOOD threshold is used for lambda2 calc (adding score to lambda2 is done before calling this function).
    // However, regions defined by minimum coverage of such good reads is not used later.
    // Because of this reason, reads having a higher score than MEDIUM are discarded here.
    // regions defined by minimum coverage of such medium good reads are not used either, but keep just in case.
    for (j = 0; j < vc.n; ++j) {
        uint32_t old_med_cov  = med_cov;
        if(vc.a[j]&2){
            if (vc.a[j]&1){
                if(vc.a[j]&4)
                    med_cov -= min_cov;
                else
                    --med_cov;
            } else {
                if(vc.a[j]&4)
                    med_cov += min_cov;
                else
                    ++med_cov;
            }
        }
        if (old_med_cov < min_cov && med_cov >= min_cov) {
            med_start = vc.a[j];
        } else if(old_med_cov >= min_cov && med_cov < min_cov) {
            uint32_t mlen = (vc.a[j]>>3) - med_start;
            if(mlen > 0){
                lq_subcoords_t mc;
                mc.start = med_start; mc.end = vc.a[j];
                kv_push(lq_subcoords_t, km, mcoords, mc);
                //printf("[DEBUG] start and end: %d, %d\n", med_start>>3, mc.end >> 3);
                lq_subcoords_t marker;
                marker.start = med_start | 0x4; //if 0x4th bit is on, cov += min_cov.
                marker.end   = vc.a[j] | 0x4; //if 0x4th bit and 0x1 are on, cov -= min_cov.
                kv_push(lq_subcoords_t, 0, *v, marker);
            }
        }
    }
    kfree(km, vc.a);
    
    int rm_tot = 0;
    for(i=0; i<cv->n; ++i){
        int flag = 0;
        if(cv->a[i].start&4){
            flag = 0; // marker
        } else {
            for(j=0; j<mcoords.n; ++j){
                if(cv->a[i].start >= mcoords.a[j].start && cv->a[i].end <= mcoords.a[j].end)
                    flag |= 1;
            }
        }
        
        if(!flag){
            kv_push(lq_subcoords_t, 0, *v, cv->a[i]);
            //printf("[DEBUG] unfiltered... start and end: %d, %d\n", cv->a[i].start>>3, cv->a[i].end >> 3);
        } else
            rm_tot++;
    }
    //printf("[DEBUG] before_filtering: %zu, after_filetering: %zu\n", cv->n, v->n);
    
    kfree(km, cv->a);
    kfree(km, mcoords.a);
}

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

int lq_map_frag_mod(const mm_idx_t *mi, int n_segs, const int *qlens, const char **seqs, const char **quals, int *n_regs, mm_reg1_t **regs,
                mm_tbuf_t *b, const mm_mapopt_t *opt, const lq_fltopt_t *fopt, int min_score_med_good, const char *qname,
                uint64_t *lambdas, uint64_t *lambdas2, lq_minimizer_cnt_v *m_cnts, lq_subcoords_v **ovlp_cors, float *ks, int off_q)
{
    int i, rep_len, qlen_sum, n_regs0, n_mini_pos;
    
    int max_chain_gap_qry, max_chain_gap_ref, is_splice = !!(opt->flag & MM_F_SPLICE), is_sr = !!(opt->flag & MM_F_SR);
    uint32_t hash;
    int64_t n_a;
    uint64_t *u, *mini_pos;
    
    m_array mc = m_cnts->a[off_q];
    
    mm128_t *a;
    mm128_v mv = {0,0,0};
    lq_subcoords_v cv = {0,0,0};
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
    
    //a = lq_chain_dp_novec(mi, max_chain_gap_ref, max_chain_gap_qry, fopt, opt->bw, opt->max_chain_skip, opt->min_cnt, opt->min_chain_score, min_score_med_good, is_splice, n_segs, n_a, a, &n_regs0, &u, b->km, qlen_sum);
    a = mm_chain_dp(max_chain_gap_ref, max_chain_gap_qry, opt->bw, opt->max_chain_skip, opt->min_cnt, opt->min_chain_score, is_splice, n_segs, n_a, a, &n_regs0, &u, b->km);
    
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
            //a = lq_chain_dp_novec(mi, max_chain_gap_ref, max_chain_gap_qry, fopt, opt->bw, opt->max_chain_skip, opt->min_cnt, opt->min_chain_score, min_score_med_good, is_splice, n_segs, n_a, a, &n_regs0, &u, b->km, qlen_sum);
            a = mm_chain_dp(max_chain_gap_ref, max_chain_gap_qry, opt->bw, opt->max_chain_skip, opt->min_cnt, opt->min_chain_score, is_splice, n_segs, n_a, a, &n_regs0, &u, b->km);
        }
    }
    
    regs0 = mm_gen_regs(b->km, hash, qlen_sum, n_regs0, u, a);
    
    chain_post(opt, max_chain_gap_ref, mi, b->km, qlen_sum, n_segs, qlens, &n_regs0, regs0, a);
    if (!is_sr) lq_cnt_match(mi, qlen_sum, n_regs0, regs0, a, n_mini_pos, mini_pos, b->km, &mc,
                             &cv, &lambdas[off_q], &lambdas2[off_q], min_score_med_good, &ks[off_q], fopt);
    //if (!is_sr) mm_est_err(mi, qlen_sum, n_regs0, regs0, a, n_mini_pos, mini_pos);
    
    //if(ovlp_cors[off_q]->n > 100)
    filter_redundant_coords(ovlp_cors[off_q], &cv, b->km, fopt->min_coverage);

    //printf("[DEBUG] cv size: %zu\n", ovlp_cors[off_q]->n);
    
    if (n_segs == 1) { // uni-segment
        regs0 = align_regs(opt, mi, b->km, qlens[0], seqs[0], quals? quals[0] : 0, &n_regs0, regs0, a);
        mm_set_mapq(n_regs0, regs0, opt->min_chain_score, opt->a, rep_len, is_sr);
        n_regs[0] = n_regs0, regs[0] = regs0;
    } else { // multi-segment
        mm_seg_t *seg;
        seg = mm_seg_gen(b->km, hash, n_segs, qlens, n_regs0, regs0, n_regs, regs, a); // split fragment chain to separate segment chains
        free(regs0);
        for (i = 0; i < n_segs; ++i) {
            mm_set_parent(b->km, opt->mask_level, n_regs[i], regs[i], opt->a * 2 + opt->b); // update mm_reg1_t::parent
            regs[i] = align_regs(opt, mi, b->km, qlens[i], seqs[i], quals? quals[i] : 0, &n_regs[i], regs[i], seg[i].a);
            mm_set_mapq(n_regs[i], regs[i], opt->min_chain_score, opt->a, rep_len, is_sr);
        }
        mm_seg_free(b->km, n_segs, seg);
        if (n_segs == 2 && opt->pe_ori >= 0 && (opt->flag&MM_F_CIGAR))
            mm_pair(b->km, max_chain_gap_ref, opt->pe_bonus, opt->a * 2 + opt->b, opt->a, qlens, n_regs, regs); // pairing
    }
    
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
        }
    }
    
    return 1;
}

/*

int lq_map_frag(const mm_idx_t *mi, int n_segs, const int *qlens, const char **seqs, const char **quals, int *n_regs, mm_reg1_t **regs,
                mm_tbuf_t *b, const mm_mapopt_t *opt, const lq_fltopt_t *fopt, int min_score_med_good, const char *qname, mm_tbuf_t *b2,
                uint64_t *lambdas, uint64_t *lambdas2, lq_subcoords_v **ovlp_cors, lq_minimizer_cnt_v **m_cnts, float *ks, int off_q)
{
    int i, rep_len, qlen_sum, n_regs0, n_mini_pos;
    
    int max_chain_gap_qry, max_chain_gap_ref, is_splice = !!(opt->flag & MM_F_SPLICE), is_sr = !!(opt->flag & MM_F_SR);
    uint32_t hash;
    int64_t n_a;
    uint64_t *u, *mini_pos;
    
    // YF memo: this CAUSED huge memory fragmentation
    if(ovlp_cors[off_q] == NULL){
        ovlp_cors[off_q] = (lq_subcoords_v *)malloc(sizeof(lq_subcoords_v));
        kv_init(*ovlp_cors[off_q]);
    }
    lq_subcoords_v *vc = ovlp_cors[off_q];
    
    // temporary holder
    uint64_t lambda_temp;
    uint64_t lambda2_temp;
    
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
        for(i=0; i<mv.n; i++){
            kv_push(uint16_t, b2->km, *m_cnts[off_q], 0);
        }
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
    a = lq_chain_dp(mi, max_chain_gap_ref, max_chain_gap_qry, fopt, opt->bw, opt->max_chain_skip, opt->min_cnt, opt->min_chain_score, min_score_med_good, is_splice, n_segs, n_a, a, &n_regs0, &u, b->km, qlen_sum, &lambdas[off_q], &lambdas2[off_q], vc);
    //a = mm_chain_dp(max_chain_gap_ref, max_chain_gap_qry, opt->bw, opt->max_chain_skip, opt->min_cnt, opt->min_chain_score, is_splice, n_segs, n_a, a, &n_regs0, &u, b->km);
    
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
            kv_init(*vc);
            a = lq_chain_dp(mi, max_chain_gap_ref, max_chain_gap_qry, fopt, opt->bw, opt->max_chain_skip, opt->min_cnt, opt->min_chain_score, min_score_med_good, is_splice, n_segs, n_a, a, &n_regs0, &u, b->km, qlen_sum, &lambdas[off_q], &lambdas2[off_q], vc);
            //a = mm_chain_dp(max_chain_gap_ref, max_chain_gap_qry, opt->bw, opt->max_chain_skip, opt->min_cnt, opt->min_chain_score, is_splice, n_segs, n_a, a, &n_regs0, &u, b->km);
        }
    }
    
    // followed relevant code from miniasm
    // the rest part is written in compute_reliable_region()
    radix_sort_32(vc->a, vc->a + vc->n);
    
    //lq_cnt_exact_match(min_score_med_good, qlen_sum, n_regs0, u, a, n_mini_pos, mini_pos, mc);
    
    regs0 = mm_gen_regs(b->km, hash, qlen_sum, n_regs0, u, a);
    
    chain_post(opt, max_chain_gap_ref, mi, b->km, qlen_sum, n_segs, qlens, &n_regs0, regs0, a);
    //if (!is_sr) lq_cnt_match(mi, qlen_sum, lambdas[off_q], n_regs0, regs0, a, n_mini_pos, mini_pos, mc, min_score_med_good, &ks[off_q], fopt);
    if (!is_sr) mm_est_err(mi, qlen_sum, n_regs0, regs0, a, n_mini_pos, mini_pos);
    
     //if (n_segs == 1) { // uni-segment
        //regs0 = align_regs(opt, mi, b->km, qlens[0], seqs[0], quals? quals[0] : 0, &n_regs0, regs0, a);
        //mm_set_mapq(n_regs0, regs0, opt->min_chain_score, opt->a, rep_len, is_sr);
        //n_regs[0] = n_regs0, regs[0] = regs0;
     //} else { // multi-segment
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
        }
    }
    
    return 1;
}

// Multi-threaded mapping

typedef struct {
    int mini_batch_size, n_processed, n_threads, n_fp;
    const mm_mapopt_t *opt;
    mm_bseq_file_t **fp;
    const mm_idx_t *mi;
    kstring_t str;
    const lq_fltopt_t *fopt; // YF memo:
    int *off_q; // YF memo:
    uint32_t min_score_med_good; //YF memo:
    uint64_t *lambdas; //YF memo:
    uint64_t *lambdas2; //YF memo:
    //float *avg_ks; //YF memo:
    //lq_subcoords_v **ovlp_cors; //YF memo:
    //lq_minimizer_cnt_v **m_cnts; // YF memo:
} pipeline_t;

typedef struct {
    const pipeline_t *p;
    int n_seq, n_frag;
    int off_th_q; // YF memo:
    mm_bseq1_t *seq;
    int *n_reg, *seg_off, *n_seg;
    mm_reg1_t **reg;
    mm_tbuf_t **buf;
    mm_tbuf_t **buf2; //YF memo: buffer for vars will be used in pipeline step == 2.
    lq_minimizer_cnt_v **m_th_cnts; // YF memo: this vector will be stored on buf2(km).
    int *seq2tid; // YF memo: need for freeing up vars in buf2
} step_t; //storages for each thred

// worker fucntion for each thread
static void worker_for(void *_data, long i, int tid) // kt_for() callback
{
    step_t *s = (step_t*)_data;
    int qlens[MM_MAX_SEG], j, off = s->seg_off[i], pe_ori = s->p->opt->pe_ori, is_sr = !!(s->p->opt->flag & MM_F_SR);
    //int off_q = off + s->off_th_q; //YF memo: off is the index within a minibatch and s->off_th_q is the offset for the chunk.
    const char *qseqs[MM_MAX_SEG], *quals[MM_MAX_SEG];
    mm_tbuf_t *b  = s->buf[tid];
    mm_tbuf_t *b2 = s->buf2[tid];
    s->seq2tid[off] = tid;
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
            mm_map_frag(s->p->mi, 1, &qlens[j], &qseqs[j], &quals[j], &s->n_reg[off+j], &s->reg[off+j], b, s->p->opt, s->seq[off+j].name);
            //lq_map_frag(s->p->mi, 1, &qlens[j], &qseqs[j], &quals[j], &s->n_reg[off+j], &s->reg[off+j], b, s->p->opt, s->p->fopt, s->p->min_score_med_good, s->seq[off+j].name, b2, s->p->lambdas, s->p->lambdas2,s->p->ovlp_cors, s->m_th_cnts, s->p->avg_ks, off);
    } else {
        mm_map_frag(s->p->mi, s->n_seg[i], qlens, qseqs, quals, &s->n_reg[off], &s->reg[off], b, s->p->opt, s->seq[off].name);
        //lq_map_frag(s->p->mi, s->n_seg[i], qlens, qseqs, quals, &s->n_reg[off], &s->reg[off], b, s->p->opt, s->p->fopt, s->p->min_score_med_good, s->seq[off].name, b2, s->p->lambdas, s->p->lambdas2, s->p->ovlp_cors, s->m_th_cnts, s->p->avg_ks, off);
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
            s->buf  = (mm_tbuf_t**)calloc(p->n_threads, sizeof(mm_tbuf_t*));
            s->buf2 = (mm_tbuf_t**)calloc(p->n_threads, sizeof(mm_tbuf_t*));
            s->seq2tid = (int*)calloc(s->n_seq, sizeof(int));
            for (i = 0; i < p->n_threads; ++i){
                s->buf[i] = mm_tbuf_init();
                s->buf2[i] = mm_tbuf_init();
            }
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
            s->m_th_cnts = (lq_minimizer_cnt_v **)calloc(s->n_seq, sizeof(lq_minimizer_cnt_v));
            s->off_th_q  = *p->off_q;
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
        // YF: add counts in this batch
        
        for (i=0; i < s->n_seq; ++i){
            //fprintf(stderr, "DEBUG: seq=%d, tid=%d, index=%d\n", i, s->seq2tid[i], *p->off_q + i);
            //if(p->m_cnts[*p->off_q + i]==NULL){
            //    p->m_cnts[*p->off_q + i] = (lq_minimizer_cnt_v *)malloc(sizeof(lq_minimizer_cnt_v));
            //    kv_init(*p->m_cnts[*p->off_q + i]);
            //    for(j=0; j<s->m_th_cnts[i]->n; ++j)
            //    kv_push(uint16_t, 0, *p->m_cnts[*p->off_q + i], 0);
            //}
            //for(j=0; j<s->m_th_cnts[i]->n;++j)
            //    p->m_cnts[*p->off_q + i]->a[j] += s->m_th_cnts[i]->a[j];
            
            if(s->buf2[s->seq2tid[i]]->km) {
                //kfree(s->buf2[s->seq2tid[i]]->km, s->m_th_cnts[i]->a);
            }
        }
        
        for (i = 0; i < p->n_threads; ++i){
            if(s->buf2[i]->km) {
                km_stat_t kmst;
                km_stat(s->buf2[i]->km, &kmst);
                assert(kmst.n_blocks == kmst.n_cores); // otherwise, there is a memory leak
                // largest flag & destroy step have to be omit; otherwise, some count information will be lost.
            }
        }
        
        for (i = 0; i < p->n_threads; ++i){
            mm_tbuf_destroy(s->buf[i]);
            mm_tbuf_destroy(s->buf2[i]);
        }
        free(s->buf); free(s->buf2);
        *(p->off_q) += s->n_seq; // add offset of query.
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
        free(s->m_th_cnts); //YF: array of pointers is freed here.
        fprintf(stderr, "[M::%s::%.3f*%.2f] mapped %d sequences. (Peak RSS: %.3f GB) \n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), s->n_seq, peakrss() / 1024.0 / 1024.0 / 1024.0);
        free(s);
    }
    return 0;
}

int lq_map_file_frag(const mm_idx_t *idx, int n_segs, const char **fn, const mm_mapopt_t *opt, const lq_fltopt_t *fopt, int min_score_med, int min_score_good, int n_threads, uint64_t *lambdas, uint64_t *lambdas2, lq_subcoords_v **ovlp_cors, lq_minimizer_cnt_v **m_cnts, float *ks, int *off)
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
    pl.lambdas = lambdas; // YF memo:
    pl.lambdas2 = lambdas2; // YF memo: for primary chain
    //pl.ovlp_cors = ovlp_cors; // YF memo:
    //pl.m_cnts = m_cnts; // YF:
    //pl.avg_ks = ks; // YF
    pl.min_score_med_good = (min_score_med << 16) | min_score_good; // YF memo
    pl.fopt = fopt;
    kt_pipeline(pl_threads, worker_pipeline, &pl, 3);
    free(pl.str.s);
    for (i = 0; i < n_segs; ++i)
        mm_bseq_close(pl.fp[i]);
    free(pl.fp);
    return 0;
}
*/

typedef struct {
    int mini_batch_size, n_processed, n_threads, n_fp;
    const mm_mapopt_t *opt;
    mm_bseq_file_t **fp;
    const mm_idx_t *mi;
    kstring_t str;
    const lq_fltopt_t *fopt; // YF memo:
    int *off_q; // YF memo:
    uint32_t min_score_med_good; //YF memo:
    uint64_t *lambdas; //YF memo:
    uint64_t *lambdas2; //YF memo:
    float *avg_ks; //YF memo:
    lq_subcoords_v **ovlp_cors; //YF memo:
    lq_minimizer_cnt_v *m_cnts; // YF memo:
} pipeline_t;

typedef struct {
    const pipeline_t *p;
    int n_seq, n_frag;
    int off_th_q; // YF memo:
    mm_bseq1_t *seq;
    int *n_reg, *seg_off, *n_seg;
    mm_reg1_t **reg;
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
            lq_map_frag_mod(s->p->mi, 1, &qlens[j], &qseqs[j], &quals[j], &s->n_reg[off+j], &s->reg[off+j], b, s->p->opt, s->p->fopt, s->p->min_score_med_good, s->seq[off+j].name, s->p->lambdas, s->p->lambdas2, s->p->m_cnts, s->p->ovlp_cors, s->p->avg_ks, off);
            //mm_map_frag(s->p->mi, 1, &qlens[j], &qseqs[j], &quals[j], &s->n_reg[off+j], &s->reg[off+j], b, s->p->opt, s->seq[off+j].name);
    } else {
        lq_map_frag_mod(s->p->mi, s->n_seg[i], qlens, qseqs, quals, &s->n_reg[off], &s->reg[off], b, s->p->opt, s->p->fopt, s->p->min_score_med_good, s->seq[off].name, s->p->lambdas, s->p->lambdas2, s->p->m_cnts, s->p->ovlp_cors, s->p->avg_ks, off);
        //mm_map_frag(s->p->mi, s->n_seg[i], qlens, qseqs, quals, &s->n_reg[off], &s->reg[off], b, s->p->opt, s->seq[off].name);
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
            s->buf  = (mm_tbuf_t**)calloc(p->n_threads, sizeof(mm_tbuf_t*));
            for (i = 0; i < p->n_threads; ++i){
                s->buf[i] = mm_tbuf_init();
            }
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

int lq_map_file_frag(const mm_idx_t *idx, int n_segs, const char **fn, const mm_mapopt_t *opt, const lq_fltopt_t *fopt, int min_score_med, int min_score_good, int n_threads, uint64_t *lambdas, uint64_t *lambdas2, lq_minimizer_cnt_v *m_cnts, lq_subcoords_v **ovlp_cors, float *ks, int *off)
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


int lq_map_file(const mm_idx_t *idx, const char *fn, const mm_mapopt_t *opt, const lq_fltopt_t *fopt, int min_score_med, int min_score_good, int n_threads, uint64_t *lambdas, uint64_t *lambdas2, lq_minimizer_cnt_v *m_cnts, lq_subcoords_v **ovlp_cors, float *ks, int *off)
{
    int fnum = (fn == NULL) ? 0 : 1;
    return lq_map_file_frag(idx, fnum, &fn, opt, fopt, min_score_med, min_score_good, n_threads, lambdas, lambdas2, m_cnts, ovlp_cors, ks, off);
}
