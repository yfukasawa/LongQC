//
//  example.h
//  minimap2_xcode
//
//  Created by Yoshinori Fukasawa on 2/5/18.
//  Copyright © 2018 Yoshinori Fukasawa. All rights reserved.
//

#include "kvec.h"
#include "minimap.h"

#ifndef example_h
#define example_h

#ifdef __cplusplus
extern "C" {
#endif

#define LQ_F_AVA 0x100000 //YF memo: not skip even if strcmp(qname, tname) > 0.
#define COVT 150

typedef struct {
    int rid;
    uint32_t length;
} lq_tchain;

typedef struct {
    uint32_t start;
    uint32_t end;
} lq_subcoords_t;

typedef struct {
    int max_overhang, min_ovlp;
    double min_ratio;
} lq_fltopt_t;

typedef kvec_t(lq_tchain) lq_tchain_v;
typedef kvec_t(uint32_t) lq_subcoords_v;
typedef kvec_t(uint32_t) lq_minimizer_cnt_v;
typedef kvec_t(lq_subcoords_t) lq_region_v;
    
int lq_cnt_exact_match(int min_sc_m_g, int qlen, int n_u, uint64_t *u, mm128_t *a, int n_mini_pos, uint64_t *mini_pos, lq_minimizer_cnt_v *m_cnts);
void lq_cnt_match(const mm_idx_t *mi, int qlen, uint64_t lambda, int n_regs, mm_reg1_t *regs, const mm128_t *a, int32_t n, const uint64_t *mini_pos, lq_minimizer_cnt_v *m_cnts, int min_sc_m_g, float *_avg_k, const lq_fltopt_t *fopt);
    
#ifdef __cplusplus
}
#endif

#endif /* example_h */
