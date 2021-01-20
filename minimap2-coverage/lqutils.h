//
//  lqutils.h
//  minimap2_xcode
//
//  Created by Yoshinori Fukasawa on 2/5/18.
//  Copyright © 2018 Yoshinori Fukasawa. All rights reserved.
//

#include "minimap2-coverage.h"

#ifndef lq_utils
#define lq_utils

#define EPS 1e-12
#define STOP 1.0e-8
#define TINY 1.0e-30

double meanQ(char * qual, int length);
int getQV(char * qual, int threshold, int length);

double incbeta(double a, double b, double x);
double inc_gamma_l(double lambda, double x);
double inc_gamma_u(double lambda, double x);
void compute_reliable_region(lq_subcoords_v *v, uint32_t min_cov, lq_subcoords_v *coords, lq_subcoords_v *mcoords);

#ifdef __cplusplus
}
#endif

#endif /* lq_utils */

