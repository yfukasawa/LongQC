//
//  kmer-dump.h
//  mm2_cov_xcode
//
//  Created by Yoshinori Fukasawa on 19/10/2019.
//  Copyright © 2019 Yoshinori Fukasawa. All rights reserved.
//

#ifndef kmer_dump_h
#define kmer_dump_h

#include <stdio.h>
#include "minimap.h"
#include "minimap2-coverage.h"

int lq_collect_occ(const mm_idx_t *idx, int n_segs, const char **fn, const mm_mapopt_t *opt,
                   int n_threads, lq_cnt_v *m_cnts, int *off);

#endif /* kmer_dump_h */
