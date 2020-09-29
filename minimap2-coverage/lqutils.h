//
//  lqutils.h
//  minimap2_xcode
//
//  Created by Yoshinori Fukasawa on 2/5/18.
//  Copyright © 2018 Yoshinori Fukasawa. All rights reserved.
//

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

#ifdef __cplusplus
}
#endif

#endif /* lq_utils */

