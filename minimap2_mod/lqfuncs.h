//
//  lqfuncs.h
//  minimap2_xcode
//
//  Created by Yoshinori Fukasawa on 1/30/18.
//  Copyright © 2018 Yoshinori Fukasawa. All rights reserved.
//

#ifndef lqfuncs_h
#define lqfuncs_h

#include "kvec.h"
#include "minimap2-coverage.h"

#ifdef __cplusplus
extern "C" {
#endif

mm128_t *lq_chain_dp(const mm_idx_t *mi, int max_dist_x, int max_dist_y, const lq_fltopt_t *fopt, int bw, int max_skip, int min_cnt, int min_sc, int min_sc_m_g, int is_cdna, int n_segs, int64_t n, mm128_t *a, int *n_u_, uint64_t **_u, void *km, int qlen, uint64_t *lambda, uint64_t *lambda2, void *vec, void *vec_coords);

#ifdef __cplusplus
}
#endif

#endif /* lqfuncs_h */
