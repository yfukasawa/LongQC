#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "lqutils.h"
#include "minimap2-coverage.h"
#include "ksort.h"

#define sort_key_32(x) (x)
KRADIX_SORT_INIT(32, uint32_t, sort_key_32, 4)

// Assumes ascii33
// nanopore sometimes exceeds Q42, so let's have more.
/*
double q2p[] = {
    1.000000000000000, 0.794328234724000, 0.630957344480000, 0.501187233627000, 0.398107170553000, 0.316227766017000,
    0.251188643151000, 0.199526231497000, 0.158489319246000, 0.125892541179000, 0.100000000000000, 0.079432823472400,
    0.063095734448000, 0.050118723362700, 0.039810717055300, 0.031622776601700, 0.025118864315100, 0.019952623149700,
    0.015848931924600, 0.012589254117900, 0.010000000000000, 0.007943282347240, 0.006309573444800, 0.005011872336270,
    0.003981071705530, 0.003162277660170, 0.002511886431510, 0.001995262314970, 0.001584893192460, 0.001258925411790,
    0.001000000000000, 0.000794328234724, 0.000630957344480, 0.000501187233627, 0.000398107170553, 0.000316227766017,
    0.000251188643151, 0.000199526231497, 0.000158489319246, 0.000125892541179, 0.000100000000000, 0.000079432823472,
    0.000063095734448,
};
*/
double q2p[] = {
    1.000000000000000,0.794328234724281,0.630957344480193,0.501187233627272,0.398107170553497,0.316227766016838,
    0.251188643150958,0.199526231496888,0.158489319246111,0.125892541179417,0.100000000000000,0.079432823472428,
    0.063095734448019,0.050118723362727,0.039810717055350,0.031622776601684,0.025118864315096,0.019952623149689,
    0.015848931924611,0.012589254117942,0.010000000000000,0.007943282347243,0.006309573444802,0.005011872336273,
    0.003981071705535,0.003162277660168,0.002511886431510,0.001995262314969,0.001584893192461,0.001258925411794,
    0.001000000000000,0.000794328234724,0.000630957344480,0.000501187233627,0.000398107170554,0.000316227766017,
    0.000251188643151,0.000199526231497,0.000158489319246,0.000125892541180,0.000100000000000,0.000079432823472,
    0.000063095734448,0.000050118723363,0.000039810717055,0.000031622776602,0.000025118864315,0.000019952623150,
    0.000015848931925,0.000012589254118,0.000010000000000,0.000007943282347,0.000006309573445,0.000005011872336,
    0.000003981071706,0.000003162277660,0.000002511886432,0.000001995262315,0.000001584893193,0.000001258925412,
    0.000001000000000,0.000000794328235,0.000000630957345,0.000000501187234,0.000000398107171,0.000000316227766,
    0.000000251188643,0.000000199526232,0.000000158489319,0.000000125892541,0.000000100000000,0.000000079432824,
    0.000000063095735,0.000000050118723,0.000000039810717,0.000000031622777,0.000000025118864,0.000000019952623,
    0.000000015848932,0.000000012589254,0.000000010000000,0.000000007943282,0.000000006309574,0.000000005011872,
    0.000000003981072,0.000000003162278,0.000000002511886,0.000000001995262,0.000000001584893,0.000000001258925,
    0.000000001000000,0.000000000794328,0.000000000630957,0.000000000501187,0.000000000398107,0.000000000316228,
    0.000000000251189,0.000000000199526,0.000000000158489,0.000000000125893,0.000000000100000,0.000000000079433,
    0.000000000063096,0.000000000050119,0.000000000039811,0.000000000031623,0.000000000025119,0.000000000019953,
    0.000000000015849,0.000000000012589,0.000000000010000,0.000000000007943,0.000000000006310,0.000000000005012,
    0.000000000003981,0.000000000003162,0.000000000002512,0.000000000001995,0.000000000001585,0.000000000001259,
    0.000000000001000,0.000000000000794,0.000000000000631,0.000000000000501,0.000000000000398,0.000000000000316,
    0.000000000000251,
};

double meanQ(char * qual, int length){
    //double a;
    int i; double sum = 0.0;
    for(i=0; i<length; i++){
        sum += q2p[(int)qual[i] - 33];
    }
    return -10*log10(sum/length);
}

/*
double meanQ(char * qual, int n){
    double a;
    int i; double sum = 0.0;
    for(i=0; i<n; i++){
        // YF memo: we should use table here. It should be faster.
        sum += pow(10.0, -1 * (double)((int)qual[i] - 33)/10.0 );
    }
    return -10*log10(sum/n);
}
*/

int getQV(char * qual, int threshold, int length){
    int i, sum = 0;
    int _t = threshold + 33;
    for(i=0; i<length; i++){
        if((int)qual[i] > _t)
	    sum++;
    }
    return sum;
}

// so dirty. need more revision
void compute_reliable_region(lq_subcoords_v *v, uint32_t min_cov, lq_subcoords_v *coords, lq_subcoords_v *mcoords)
{
    int j;
    uint32_t start, cov;
    uint32_t med_start, med_cov;
    start = cov = med_start = med_cov = 0;
    
    kvec_t(uint32_t) vc = {0,0,0};
    for(j=0; j<v->n;++j){
        kv_push(uint32_t, 0, vc, v->a[j].start);
        kv_push(uint32_t, 0, vc, v->a[j].end);
    }
    radix_sort_32(vc.a, vc.a + vc.n);
    
    for (j = 0, cov = 0; j < vc.n; ++j) {
        uint32_t old_cov      = cov;
        uint32_t old_med_cov  = med_cov;
        if (vc.a[j]&1){
            --cov;
            if(vc.a[j]&2){
                if(vc.a[j]&4){
                    med_cov -= min_cov;
                    cov -= (min_cov - 1); //too much culling;
                }
                else
                    --med_cov;
            }
        }
        else {
            ++cov;
            if(vc.a[j]&2){
                if(vc.a[j]&4){
                    med_cov += min_cov;
                    cov     += (min_cov - 1); // too much
                }
                else
                    ++med_cov;
            }
        }
        if (old_cov < min_cov && cov >= min_cov) {
            start = vc.a[j]>>3;
            // just in case. cov and med_cov and good_cov is not 100% exclusive.
            if(old_med_cov < min_cov && med_cov >= min_cov)
                med_start = vc.a[j]>>3;
        } else if (old_cov >= min_cov && cov < min_cov) {
            uint32_t len = (vc.a[j]>>3) - start;
            if (len > 0){
                lq_subcoords_t c;
                c.start = start; c.end = vc.a[j]>>3;
                kv_push(lq_subcoords_t, 0, *coords, c);
            }
            // just in case. cov and med_cov and good_cov is not 100% exclusive.
            if (old_med_cov >= min_cov && med_cov < min_cov){
                uint32_t mlen = (vc.a[j]>>3) - med_start;
                if(mlen > 0){
                    lq_subcoords_t mc;
                    mc.start = med_start; mc.end = vc.a[j]>>3;
                    kv_push(lq_subcoords_t, 0, *mcoords, mc);
                }
            }
        } else if (old_med_cov < min_cov && med_cov >= min_cov) {
            med_start = vc.a[j]>>3;
        } else if(old_med_cov >= min_cov && med_cov < min_cov) {
            uint32_t mlen = (vc.a[j]>>3) - med_start;
            if(mlen > 0){
                lq_subcoords_t mc;
                mc.start = med_start; mc.end = vc.a[j]>>3;
                kv_push(lq_subcoords_t, 0, *mcoords, mc);
            }
        }
    }
    kv_destroy(vc);
}


//// so dirty. need more revision
//void compute_reliable_region(lq_subcoords_v *vc, uint32_t min_cov, lq_region_v *coords, lq_region_v *mcoords, lq_region_v *gcoords)
//{
//    int j;
//    uint32_t start, cov;
//    uint32_t med_start, med_cov;
//    uint32_t good_start, good_cov;
//    start = cov = med_start = med_cov = good_start = good_cov = 0;
//    
//    // followed relevant code from miniasm
//    // the rest part is written in compute_reliable_region()
//    radix_sort_32(vc->a, vc->a + vc->n);
//    
//    for (j = 0, cov = 0; j < vc->n; ++j) {
//        uint32_t old_cov      = cov;
//        uint32_t old_med_cov  = med_cov;
//        uint32_t old_good_cov = good_cov;
//        if (vc->a[j]&1){
//            --cov;
//            if(vc->a[j]&2)
//                --med_cov;
//            if(vc->a[j]&4)
//                --good_cov;
//        }
//        else {
//            ++cov;
//            if(vc->a[j]&2)
//                ++med_cov;
//            if(vc->a[j]&4)
//                ++good_cov;
//        }
//        if (old_cov < min_cov && cov >= min_cov) {
//            start = vc->a[j]>>3;
//            // just in case. cov and med_cov and good_cov is not 100% exclusive.
//            if(old_med_cov < min_cov && med_cov >= min_cov){
//                med_start = vc->a[j]>>3;
//                if (old_good_cov < min_cov && good_cov >= min_cov)
//                    good_start = vc->a[j]>>3;
//            }
//        } else if (old_cov >= min_cov && cov < min_cov) {
//            uint32_t len = (vc->a[j]>>3) - start;
//            if (len > 0){
//                lq_subcoords_t c;
//                c.start = start; c.end = vc->a[j]>>3;
//                kv_push(lq_subcoords_t, 0, *coords, c);
//            }
//            // just in case. cov and med_cov and good_cov is not 100% exclusive.
//            if (old_med_cov >= min_cov && med_cov < min_cov){
//                uint32_t mlen = (vc->a[j]>>3) - med_start;
//                if(mlen > 0){
//                    lq_subcoords_t mc;
//                    mc.start = med_start; mc.end = vc->a[j]>>3;
//                    kv_push(lq_subcoords_t, 0, *mcoords, mc);
//                }
//                if(old_good_cov >= min_cov && good_cov < min_cov){
//                    uint32_t glen = (vc->a[j]>>3) - good_start;
//                    if(glen > 0){
//                        lq_subcoords_t gc;
//                        gc.start = good_start; gc.end = vc->a[j]>>3;
//                        kv_push(lq_subcoords_t, 0, *gcoords, gc);
//                    }
//                }
//            }
//        } else if (old_med_cov < min_cov && med_cov >= min_cov) {
//            med_start = vc->a[j]>>3;
//            // just in case. med_cov and good_cov is not 100% exclusive.
//            if (old_good_cov < min_cov && good_cov >= min_cov)
//                good_start = vc->a[j]>>3;
//        } else if(old_med_cov >= min_cov && med_cov < min_cov) {
//            uint32_t mlen = (vc->a[j]>>3) - med_start;
//            if(mlen > 0){
//                lq_subcoords_t mc;
//                mc.start = med_start; mc.end = vc->a[j]>>3;
//                kv_push(lq_subcoords_t, 0, *mcoords, mc);
//            }
//            // just in case. med_cov and good_cov is not 100% exclusive.
//            if(old_good_cov >= min_cov && good_cov < min_cov){
//                uint32_t glen = (vc->a[j]>>3) - good_start;
//                if(glen > 0){
//                    lq_subcoords_t gc;
//                    gc.start = good_start; gc.end = vc->a[j]>>3;
//                    kv_push(lq_subcoords_t, 0, *gcoords, gc);
//                }
//            }
//        } else if (old_good_cov < min_cov && good_cov >= min_cov) {
//            good_start = vc->a[j]>>3;
//        } else if (old_good_cov >= min_cov && good_cov < min_cov){
//            uint32_t glen = (vc->a[j]>>3) - good_start;
//            if(glen > 0){
//                lq_subcoords_t gc;
//                gc.start = good_start; gc.end = vc->a[j]>>3;
//                kv_push(lq_subcoords_t, 0, *gcoords, gc);
//            }
//        }
//    }
//}


/* End of messy functions*/


// cdf of poission is upper regularized gamma function.
// translate from javascript code (https://github.com/compute-io/gammainc/blob/677b930c1f1b25222368009a52dd4f0a8729d4b5/lib/number.js)
double inc_gamma_u(double lambda, double x)
{
    if(lambda <= 1 || lambda <= x)
        return 1.0 - inc_gamma_l(lambda, x);
    
    double f = 1.0 + lambda - x, C = f, D = 0.0, a, b, chg;
    int i;
    
    for(i = 1; i<10000; i++){
        a = (double)i * (x - (double)i);
        b = (double)(i << 1) + 1.0 + lambda - x;
        D = b + a * D;
        C = b + a / C;
        D = 1.0 / D;
        chg = C * D;
        f *= chg;
        
        if ( fabs( chg - 1.0 ) < EPS ) {
            break;
        }
    }
    
    return exp(x * log(lambda) - lambda - lgamma(x) - log(f) );
    
}

double inc_gamma_l(double lambda, double x)
{
    if(lambda == 0)
        return 0.0;
    
    if(lambda < 0 || lambda > x)
        return 1 - inc_gamma_u(lambda, x);
    
    double ft, r = x, c = 1.0, pws = 1.0;
    ft =  x * log(lambda) - lambda - lgamma(x);
    ft = exp(ft);
    
    do {
        r   += 1.0;
        c   *= lambda/r;
        pws += c;
    }  while( c / pws > EPS );
    
    return pws * ft / x;
}

double cdf_poisson(int k, double lambda)
{
    return inc_gamma_u(lambda, k+1);
}

/*
double cdf_neg_binom(double p, double k, double r)
{
    return 1.0 - incbeta(k+1.0, r, p);
}
*/

/*
 * zlib License
 *
 * Regularized Incomplete Beta Function
 *
 * Copyright (c) 2016, 2017 Lewis Van Winkle
 * http://CodePlea.com
 *
 * This software is provided 'as-is', without any express or implied
 * warranty. In no event will the authors be held liable for any damages
 * arising from the use of this software.
 *
 * Permission is granted to anyone to use this software for any purpose,
 * including commercial applications, and to alter it and redistribute it
 * freely, subject to the following restrictions:
 *
 * 1. The origin of this software must not be misrepresented; you must not
 *    claim that you wrote the original software. If you use this software
 *    in a product, an acknowledgement in the product documentation would be
 *    appreciated but is not required.
 * 2. Altered source versions must be plainly marked as such, and must not be
 *    misrepresented as being the original software.
 * 3. This notice may not be removed or altered from any source distribution.
 */


double incbeta(double a, double b, double x) {
    if (x < 0.0 || x > 1.0) return 1.0/0.0;
    
    /*The continued fraction converges nicely for x < (a+1)/(a+b+2)*/
    if (x > (a+1.0)/(a+b+2.0)) {
        return (1.0-incbeta(b,a,1.0-x)); /*Use the fact that beta is symmetrical.*/
    }
    
    /*Find the first part before the continued fraction.*/
    const double lbeta_ab = lgamma(a)+lgamma(b)-lgamma(a+b);
    const double front = exp(log(x)*a+log(1.0-x)*b-lbeta_ab) / a;
    
    /*Use Lentz's algorithm to evaluate the continued fraction.*/
    double f = 1.0, c = 1.0, d = 0.0;
    
    int i, m;
    for (i = 0; i <= 200; ++i) {
        m = i/2;
        
        double numerator;
        if (i == 0) {
            numerator = 1.0; /*First numerator is 1.0.*/
        } else if (i % 2 == 0) {
            numerator = (m*(b-m)*x)/((a+2.0*m-1.0)*(a+2.0*m)); /*Even term.*/
        } else {
            numerator = -((a+m)*(a+b+m)*x)/((a+2.0*m)*(a+2.0*m+1)); /*Odd term.*/
        }
        
        /*Do an iteration of Lentz's algorithm.*/
        d = 1.0 + numerator * d;
        if (fabs(d) < TINY) d = TINY;
        d = 1.0 / d;
        
        c = 1.0 + numerator / c;
        if (fabs(c) < TINY) c = TINY;
        
        const double cd = c*d;
        f *= cd;
        
        /*Check for stop.*/
        if (fabs(1.0-cd) < STOP) {
            return front * (f-1.0);
        }
    }
    
    return 1.0/0.0; /*Needed more loops, did not converge.*/
}
