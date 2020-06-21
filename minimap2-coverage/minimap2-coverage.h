//
//  example.h
//  ode
//
//  Created by Yoshinori Fukasawa on 2/5/18.
//  Copyright © 2018 Yoshinori Fukasawa. All rights reserved.
//

#ifndef minimap2_coverage
#define minimap2_coverage

#ifdef __cplusplus
extern "C" {
#endif

#include "kvec.h"
#include "minimap.h"

#define LQ_F_AVA 0x100000 //YF memo: not skip even if strcmp(qname, tname) > 0.
#define COVT 150

typedef struct {
    uint32_t start;
    uint32_t end;
} lq_subcoords_t;

typedef struct {
    int max_overhang, min_ovlp, min_coverage;
    double min_ratio;
} lq_fltopt_t;
    
typedef struct {
    uint32_t n;
    uint16_t *a;
} m_array;

typedef kvec_t(lq_subcoords_t) lq_subcoords_v;
typedef kvec_t(m_array) lq_minimizer_cnt_v;
    
//int lq_cnt_exact_match(int min_sc_m_g, int qlen, int n_u, uint64_t *u, mm128_t *a, int n_mini_pos, uint64_t *mini_pos, lq_minimizer_cnt_v *m_cnts);
int lq_map_file(const mm_idx_t *idx, const char *fn, const mm_mapopt_t *opt, const lq_fltopt_t *fopt, int min_score_med, int min_score_good,
                int n_threads, uint64_t *lambdas, uint64_t *lambdas2, lq_minimizer_cnt_v *m_cnts, lq_subcoords_v **ovlp_cors, float *ks, int *off);

void lq_cnt_match(const mm_idx_t *mi, int qlen, int n_regs, mm_reg1_t *regs, const mm128_t *a, int32_t n, const uint64_t *mini_pos, void *km, 
                  m_array *m_cnts, lq_subcoords_v *ovlp_cors, uint64_t *lambda, uint64_t *lambda2, int min_sc_m_g, float *_avg_k, const lq_fltopt_t *fopt);

void radix_sort_32(uint32_t *beg, uint32_t *end);
    
#ifdef __cplusplus
}
#endif

#endif
