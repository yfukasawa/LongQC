#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "mmpriv.h"
#include "minimap2-coverage.h"

/*
 https://stackoverflow.com/questions/1340729/how-do-you-generate-a-random-double-uniformly-distributed-between-0-and-1-from-c
 let's try this simple one.
 */
double randZeroToOne()
{
    return rand() / (RAND_MAX + 1.);
}

static inline int32_t get_for_qpos(int32_t qlen, const mm128_t *a)
{
	int32_t x = (int32_t)a->y;
	int32_t q_span = a->y>>32 & 0xff;
	if (a->x>>63)
		x = qlen - 1 - (x + 1 - q_span); // revert the position to the forward strand of query
	return x;
}

static int get_mini_idx(int qlen, const mm128_t *a, int32_t n, const uint64_t *mini_pos)
{
	int32_t x, L = 0, R = n - 1;
	x = get_for_qpos(qlen, a);
	while (L <= R) { // binary search
		int32_t m = ((uint64_t)L + R) >> 1;
		int32_t y = (int32_t)mini_pos[m];
		if (y < x) L = m + 1;
		else if (y > x) R = m - 1;
		else return m;
	}
	return -1;
}

/*
 // too buggy. massive fix is required to work on.
int lq_cnt_exact_match(int min_sc_m_g, int qlen, int n_u, uint64_t *u, mm128_t *a, int n_mini_pos, uint64_t *mini_pos, lq_minimizer_cnt_v *m_cnts) // countup exact matches on a query
{
    uint16_t min_sc_m, min_sc_g;
    int32_t n_seed, st, i, j, k, offset = 0;
    min_sc_m = (uint16_t)(min_sc_m_g >> 16);
    min_sc_g = (uint16_t) min_sc_m_g;
    int th = min_sc_m, rev;
    
    if (n_u == 0) return 0;
    
    for (i = k = 0; i < n_u; ++i) {
        if((int32_t)(u[i] >> 32) < th){
            offset += (int32_t)u[i];
            continue;
        }
        rev = a[offset].x>>63;
        n_seed = (int32_t)u[i];
        st = get_mini_idx(qlen, rev? &a[offset + n_seed - 1] : &a[offset], n_mini_pos, mini_pos);
        for (k = 1, j = st + 1; j < n_mini_pos && k < n_seed; ++j) {
            int32_t x;
            x = get_for_qpos(qlen, rev? &a[offset + n_seed - 1 - k] : &a[offset + k]);
            if (x == (int32_t)mini_pos[j])
                m_cnts->a[j] += 1;
        }
        offset += (int32_t)u[i];
    }
    return 1;
}
*/

void lq_cnt_match(const mm_idx_t *mi, int qlen, int n_regs, mm_reg1_t *regs, const mm128_t *a, int32_t n, const uint64_t *mini_pos, void *km,
                  m_array *m_cnts, lq_subcoords_v *cv, uint64_t *lambda, uint64_t *lambda2, int min_sc_m_g, float *_avg_k, const lq_fltopt_t *fopt)
{
    int i;
    uint64_t sum_k = 0;
    uint32_t qs, qe, rs, re, rl, hang5, hang3;
    
    uint16_t min_sc_m;
    min_sc_m = (uint16_t)(min_sc_m_g >> 16);
    uint16_t min_sc_g;
    min_sc_g = (uint16_t) min_sc_m_g;
    double frac = 1.0;
    
    if (n == 0) return;
    
    if(*lambda/qlen > COVT && *_avg_k != 0.0){
        return; // already reached
    } else if (*lambda/qlen > COVT && *_avg_k == 0.0) {
        frac = (double) COVT/(*lambda/qlen);
    }
    
    if(*_avg_k == 0.0){
        for (i = 0; i < n; ++i)
            sum_k += mini_pos[i] >> 32 & 0xff;
        *_avg_k = (float)sum_k / n;
    }
    
    for (i = 0; i < n_regs; ++i) {
        if(randZeroToOne() >= frac) continue;
        mm_reg1_t *r = &regs[i];
        int32_t st, j, k;
        r->div = -1.0f;
        int flag = 0;
        if (r->cnt == 0) continue;
        st = get_mini_idx(qlen, r->rev? &a[r->as + r->cnt - 1] : &a[r->as], n, mini_pos);
        if (st < 0) {
            if (mm_verbose >= 2)
                fprintf(stderr, "[WARNING] logic inconsistency in mm_est_err(). Please contact the developer.\n");
            continue;
        }
        // fragment filering
        rl = mi->seq[r->rid].len;
        qs = r->qs; qe = r->qe;
        rs = r->rs; re = r->re;
        hang5 = qs < rs ? qs : rs;
        hang3 = qlen - qe < rl - re ? qlen - qe : rl - re;
        if( (qe - qs) <  (qe - qs + hang5 + hang3) * fopt->min_ratio || hang5 > fopt->max_overhang || hang3 > fopt->max_overhang )
            continue;
        lq_subcoords_t s;
        *lambda += (qe - qs + 1);
        if(r->score0 >= min_sc_m) flag = flag|0x2;
        s.start = qs<<3|flag;
        flag = flag|0x1;
        s.end = qe<<3|flag;
        kv_push(lq_subcoords_t, km, *cv, s);
        if(r->score0 < min_sc_g) continue; // lambda2 and cnts have at least min_sc_g scores.
        //flag = flag|0x4; //YF memo: good coordinates are not in use.
        *lambda2 += (qe - qs + 1);
        if(m_cnts->a[st] < UINT16_MAX) m_cnts->a[st]++;
        for (k = 1, j = st + 1; j < n && k < r->cnt; ++j) {
            int32_t x;
            x = get_for_qpos(qlen, r->rev? &a[r->as + r->cnt - 1 - k] : &a[r->as + k]);
            if (x == (int32_t)mini_pos[j]){
                ++k;
                if(m_cnts->a[st] < UINT16_MAX) m_cnts->a[j]++;
            }
        }
    }
}

void mm_est_err(const mm_idx_t *mi, int qlen, int n_regs, mm_reg1_t *regs, const mm128_t *a, int32_t n, const uint64_t *mini_pos)
{
	int i;
	uint64_t sum_k = 0;
	float avg_k;

	if (n == 0) return;
	for (i = 0; i < n; ++i)
		sum_k += mini_pos[i] >> 32 & 0xff;
	avg_k = (float)sum_k / n;

	for (i = 0; i < n_regs; ++i) {
		mm_reg1_t *r = &regs[i];
		int32_t st, en, j, k, n_match, n_tot, l_ref;
		r->div = -1.0f;
		if (r->cnt == 0) continue;
		st = en = get_mini_idx(qlen, r->rev? &a[r->as + r->cnt - 1] : &a[r->as], n, mini_pos);
		if (st < 0) {
			if (mm_verbose >= 2)
				fprintf(stderr, "[WARNING] logic inconsistency in mm_est_err(). Please contact the developer.\n");
			continue;
		}
		l_ref = mi->seq[r->rid].len;
		for (k = 1, j = st + 1, n_match = 1; j < n && k < r->cnt; ++j) {
			int32_t x;
			x = get_for_qpos(qlen, r->rev? &a[r->as + r->cnt - 1 - k] : &a[r->as + k]);
			if (x == (int32_t)mini_pos[j])
                ++k, en = j, ++n_match;
		}
		n_tot = en - st + 1;
		if (r->qs > avg_k && r->rs > avg_k) ++n_tot;
		if (qlen - r->qs > avg_k && l_ref - r->re > avg_k) ++n_tot;
        //fprintf(stderr, "[DEBUG] n_match %d\n", n_match);
		r->div = logf((float)n_tot / n_match) / avg_k;
	}
}
