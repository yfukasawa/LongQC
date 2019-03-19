#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include "minimap.h"
#include "mmpriv.h"
#include "kalloc.h"
#include "kvec.h"
#include "minimap2-coverage.h"

static const char LogTable256[256] = {
#define LT(n) n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n
    -1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
    LT(4), LT(5), LT(5), LT(6), LT(6), LT(6), LT(6),
    LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7)
};

static inline int ilog2_32(uint32_t v)
{
    register uint32_t t, tt;
    if ((tt = v>>16)) return (t = tt>>8) ? 24 + LogTable256[t] : 16 + LogTable256[tt];
    return (t = v>>8) ? 8 + LogTable256[t] : LogTable256[v];
}

/*
 int32_t i:   //first seed index in the chain
 int32_t j:   //last seed index in the chain
 int32_t ql:  //query length
 void *vec:   //pointer to struct kvec_t(lq_tchain)
 */
void lq_compute_fuzzy_overlap(const mm_idx_t *mi, const lq_fltopt_t *fopt, mm128_t *a, int32_t i, int32_t j, int is_bet_med, int is_bet_good, int32_t ql, lq_tchain_v *v, lq_subcoords_v *vc, uint64_t *lambda, uint64_t *lambda2)
{
    int flag = 0;
    int32_t m_span;
    int32_t qs, qe, rid, rs, re, rl, hang5, hang3;
    lq_tchain *c;
    kv_pushp(lq_tchain, 0, *v, &c);
    
    /*
    int max_overhang = 2000;
    int min_ovlp     = 1000;
    double  internal_ratio = 0.4;
     */
    
    rid    = a[i].x<<1>>33;
    rl     = mi->seq[rid].len;
    m_span = (int32_t)(a[i].y>>32&0xff);
    
    rs = (int32_t)a[j].x + 1 > m_span? (int32_t)a[j].x + 1 - m_span : 0; // copied from set_coord
    re = (int32_t)a[i].x + 1;
    
    // in case of rev comp match, 1st bit of a[i].x should be flagged
    if(a[i].x>>63){
        qs = ql - (a[i].y + 1);
        qe = ql - (a[j].y + 1 - m_span);
    } else {
        qs = a[j].y + 1 - m_span;
        qe = a[i].y + 1;
    }
    
    if(qs < rs){
        hang5 = qs;
    } else {
        hang5 = rs;
    }
    
    if(ql - qe < rl - re){
        hang3 = ql - qe;
    } else {
        hang3 = rl - re;
    }
    
    //fprintf(stderr, "DEBUG\tquery_pos\t%d\t%d\ttarget_pos\t%d\t%d\t%d\thang5,3\t%d\t%d\n", qs, qe, mi->seq[rid].len, rs, re, hang5, hang3);

    if( (qe - qs) < fopt->min_ovlp || (re - rs) < fopt->min_ovlp ){
        return;
    }

    if( (qe - qs) <  (qe - qs + hang5 + hang3) * fopt->min_ratio || hang5 > fopt->max_overhang || hang3 > fopt->max_overhang ){
        return;
    }
    
    c->rid    = rid;
    c->length = re - rs + 1;
    
    // YF memo: follow the way of miniasm.
    *lambda += (qe - qs + 1);
    if(is_bet_med)
        flag = flag|0x2;
    if(is_bet_good){
        flag = flag|0x4;
        *lambda2 += (qe - qs + 1);
    }
    kv_push(uint32_t, 0, *vc, qs<<3|flag);
    flag = flag|0x1;
    kv_push(uint32_t, 0, *vc, qe<<3|flag);
}

mm128_t *lq_chain_dp(const mm_idx_t *mi, int max_dist_x, int max_dist_y, const lq_fltopt_t *fopt, int bw, int max_skip, int min_cnt, int min_sc, int min_sc_m_g, int is_cdna, int n_segs, int64_t n, mm128_t *a, int *n_u_, uint64_t **_u, void *km, int qlen, uint64_t *lambda, uint64_t *lambda2, lq_tchain_v *_vec, lq_subcoords_v *_vec_coords)
{ // TODO: make sure this works when n has more than 32 bits
    
    uint16_t min_sc_m, min_sc_g;
    int32_t k, *f, *p, *t, *v, n_u, n_v;
    int64_t i, j, st = 0, l;
    uint64_t *u, *u2, sum_qspan = 0;
    float avg_qspan;
    mm128_t *b, *w;
    
    min_sc_m = (uint16_t)(min_sc_m_g >> 16);
    min_sc_g = (uint16_t) min_sc_m_g;
    
    if (_u) *_u = 0, *n_u_ = 0;
    f = (int32_t*)kmalloc(km, n * 4);
    p = (int32_t*)kmalloc(km, n * 4);
    t = (int32_t*)kmalloc(km, n * 4);
    v = (int32_t*)kmalloc(km, n * 4);
    memset(t, 0, n * 4);
    
    for (i = 0; i < n; ++i) sum_qspan += a[i].y>>32&0xff;
    avg_qspan = (float)sum_qspan / n;
    
    // fill the score and backtrack arrays
    for (i = 0; i < n; ++i) {
        uint64_t ri = a[i].x;
        int64_t max_j = -1;
        int32_t qi = (int32_t)a[i].y, q_span = a[i].y>>32&0xff; // NB: only 8 bits of span is used!!!
        int32_t max_f = q_span, n_skip = 0, min_d;
        int32_t sidi = (a[i].y & MM_SEED_SEG_MASK) >> MM_SEED_SEG_SHIFT;
        while (st < i && ri - a[st].x > max_dist_x) ++st;
        for (j = i - 1; j >= st; --j) {
            int64_t dr = ri - a[j].x;
            int32_t dq = qi - (int32_t)a[j].y, dd, sc, log_dd;
            int32_t sidj = (a[j].y & MM_SEED_SEG_MASK) >> MM_SEED_SEG_SHIFT;
            if ((sidi == sidj && dr == 0) || dq <= 0) continue; // don't skip if an anchor is used by multiple segments; see below
            if ((sidi == sidj && dq > max_dist_y) || dq > max_dist_x) continue;
            dd = dr > dq? dr - dq : dq - dr;
            if (sidi == sidj && dd > bw) continue;
            if (n_segs > 1 && !is_cdna && sidi == sidj && dr > max_dist_y) continue;
            min_d = dq < dr? dq : dr;
            sc = min_d > q_span? q_span : dq < dr? dq : dr;
            log_dd = dd? ilog2_32(dd) : 0;
            if (is_cdna || sidi != sidj) {
                int c_log, c_lin;
                c_lin = (int)(dd * .01 * avg_qspan);
                c_log = log_dd;
                if (sidi != sidj && dr == 0) ++sc; // possibly due to overlapping paired ends; give a minor bonus
                else if (dr > dq || sidi != sidj) sc -= c_lin < c_log? c_lin : c_log;
                else sc -= c_lin + (c_log>>1);
            } else sc -= (int)(dd * .01 * avg_qspan) + (log_dd>>1);
            sc += f[j];
            if (sc > max_f) {
                max_f = sc, max_j = j;
                if (n_skip > 0) --n_skip;
            } else if (t[j] == i) {
                if (++n_skip > max_skip)
                    break;
            }
            if (p[j] >= 0) t[p[j]] = i;
        }
        f[i] = max_f, p[i] = max_j;
        v[i] = max_j >= 0 && v[max_j] > max_f? v[max_j] : max_f; // v[] keeps the peak score up to i; f[] is the score ending at i, not always the peak
    }
    
    // find the ending positions of chains
    memset(t, 0, n * 4);
    for (i = 0; i < n; ++i)
        if (p[i] >= 0) t[p[i]] = 1;
    for (i = n_u = 0; i < n; ++i)
        if (t[i] == 0 && v[i] >= min_sc)
            ++n_u;
    if (n_u == 0) {
        kfree(km, a); kfree(km, f); kfree(km, p); kfree(km, t); kfree(km, v);
        return 0;
    }
    
    u = (uint64_t*)kmalloc(km, n_u * 8);
    for (i = n_u = 0; i < n; ++i) {
        if (t[i] == 0 && v[i] >= min_sc) {
            j = l = i;
            while (j >= 0 && f[j] < v[j]) j = p[j]; // find the peak that maximizes f[]
            if (j < 0){
                j = i; // TODO: this should really be assert(j>=0)
            } else if (i != j){
                //printf("debug: j >= 0 in chain_lq %llu %llu\n", j, i);
            }
            u[n_u++] = (uint64_t)f[j] << 32 | j; // YF memo: f[] is the chain score, and j is basically i
            
        }
    }
    
    radix_sort_64(u, u + n_u);
    for (i = 0; i < n_u>>1; ++i) { // reverse, s.t. the highest scoring chain is the first
        uint64_t t = u[i];
        u[i] = u[n_u - i - 1], u[n_u - i - 1] = t;
    }
    
    // backtrack
    memset(t, 0, n * 4);
    for (i = n_v = k = 0; i < n_u; ++i) { // starting from the highest score
        int32_t n_v0 = n_v, k0 = k; int64_t j_prev = -1; int64_t seed_st; int sc=0;
        j = seed_st = (int32_t)u[i];
        do {
            v[n_v++] = j;
            t[j] = 1;
            j_prev = j;
            j = p[j];
        } while (j >= 0 && t[j] == 0);
        
        if (j < 0) {
            if (n_v - n_v0 >= min_cnt){
                u[k++] = u[i]>>32<<32 | (n_v - n_v0);
                lq_compute_fuzzy_overlap(mi, fopt, a, seed_st, j_prev, u[i]>>32 >= min_sc_m, u[i]>>32 >= min_sc_g, qlen, _vec, _vec_coords, lambda, lambda2);
            }
        } else if ((int32_t)(u[i]>>32) - f[j] >= min_sc) {
            if (n_v - n_v0 >= min_cnt){
                sc = (int32_t)(u[i]>>32) - f[j];
                u[k++] = (int64_t) sc << 32 | (n_v - n_v0);
                lq_compute_fuzzy_overlap(mi, fopt, a, seed_st, j_prev, sc >= min_sc_m, sc >= min_sc_g, qlen, _vec, _vec_coords, lambda, lambda2);
            }
        }
        if (k0 == k) n_v = n_v0; // no new chain added, reset
    }
    *n_u_ = n_u = k, *_u = u; // NB: note that u[] may not be sorted by score here
    
    // free temporary arrays
    kfree(km, f); kfree(km, p); kfree(km, t);
    
    // write the result to b[]
    b = (mm128_t*)kmalloc(km, n_v * sizeof(mm128_t));
    for (i = 0, k = 0; i < n_u; ++i) {
        int32_t k0 = k, ni = (int32_t)u[i];
        for (j = 0; j < ni; ++j)
            b[k] = a[v[k0 + (ni - j - 1)]], ++k;
    }
    kfree(km, v);
    
    // sort u[] and a[] by a[].x, such that adjacent chains may be joined (required by mm_join_long)
    w = (mm128_t*)kmalloc(km, n_u * sizeof(mm128_t));
    for (i = k = 0; i < n_u; ++i) {
        w[i].x = b[k].x, w[i].y = (uint64_t)k<<32|i;
        k += (int32_t)u[i];
    }
    radix_sort_128x(w, w + n_u);
    u2 = (uint64_t*)kmalloc(km, n_u * 8);
    for (i = k = 0; i < n_u; ++i) {
        int32_t j = (int32_t)w[i].y, n = (int32_t)u[j];
        u2[i] = u[j];
        memcpy(&a[k], &b[w[i].y>>32], n * sizeof(mm128_t));
        k += n;
    }
    memcpy(u, u2, n_u * 8);
    memcpy(b, a, k * sizeof(mm128_t)); // write _a_ to _b_ and deallocate _a_ because _a_ is oversized, sometimes a lot
    kfree(km, a); kfree(km, w); kfree(km, u2);
    return b;
}
