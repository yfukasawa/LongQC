import os, re, logging
#import bindings.python.edlib
import numpy as np
import edlib

def _cutr(reads, adp, th, r, *, len_list=[]):
    iden_max  = -1
    match_num = 0
    cut_pos   = []
    repat = re.compile(r'(\d+)[DHIMNPSX=]{1}')
    if len(reads[0]) > 2:
        has_qual = True
    else:
        has_qual = False
    for read in reads:
        if not len_list:
            pass
        else:
            len_list.append( len(read[1]) )
        result = edlib.align(adp, read[1][-r:], mode="HW", task='path')
        identity = 1.0 - float(result['editDistance']/np.sum([int(i) for i in repat.findall(result['cigar'])]))
        if identity >th:
            #print(read[0], result, identity, result['locations'], result['cigar'], identity)
            start = len(read[1]) - r + result['locations'][0][0]
            cut_pos.append(r-result['locations'][0][0])
            match_num += 1        
            if identity > iden_max:
                iden_max = identity
            read[1] = read[1][:start]
            if has_qual:
                read[2] = read[2][:start]

    return (iden_max, match_num, cut_pos)

def _cutf(reads, adp, th, r, *, len_list=[]):
    iden_max  = -1
    match_num = 0
    cut_pos   = []
    repat = re.compile(r'(\d+)[DHIMNPSX=]{1}')
    if len(reads[0]) > 2:
        has_qual = True
    else:
        has_qual = False
    for read in reads:
        if not len_list:
            pass
        else:
            len_list.append( len(read[1]) )
        result = edlib.align(adp, read[1][:r], mode="HW", task='path')
        identity = 1.0 - float(result['editDistance']/np.sum([int(i) for i in repat.findall(result['cigar'])]))
        if identity >th:
            #print(read[0], result, identity, result['locations'], result['cigar'], identity)
            end = result['locations'][0][1]
            cut_pos.append(end)
            match_num += 1        
            if identity > iden_max:
                iden_max = identity
            read[1] = read[1][end+1:]
            if has_qual:
                read[2] = read[2][end+1:]

    return (iden_max, match_num, cut_pos)

def cut_adapter(reads, *, len_list=None, adp_t=None, adp_b = None, logger=None, th=0.75, length=150):
    iden_max5 = 0.0
    iden_max3 = 0.0

    if not adp_t and not adp_b:
        logger.WARNING()
        return None

    if adp_t and adp_b:
        (iden_max5, mnum5, cpos5) = _cutf(reads, adp_t, th, length, len_list=len_list)
        (iden_max3, mnum3, cpos3) = _cutr(reads, adp_b, th, length)
        logger.info("Adapter Sequence: %s, max identity:%f and the number of trimmed reads: %d" % (adp_t, iden_max5, mnum5))
        logger.info("Adapter Sequence: %s, max identity:%f and the number of trimmed reads: %d" % (adp_b, iden_max3, mnum3))
        return ((iden_max5, mnum5, cpos5), (iden_max3, mnum3, cpos3))

    if adp_t and not adp_b:
        (iden_max5, mnum5, cpos5) = _cutf(reads, adp_t, th, length, len_list=len_list)
        logger.info("Adapter Sequence: %s, max identity:%f and the number of trimmed reads: %d" % (adp_t, iden_max5, mnum5))
        return (iden_max5, mnum, cpos)

    if not adp_t and adp_b:
        (iden_max3, mnum3, cpos3) = _cutr(reads, adp_b, th, length, len_list=len_list)
        logger.info("Adapter Sequence: %s, max identity:%f and the number of trimmed reads: %d" % (adp_b, iden_max3, mnum3))
        return (iden_max3, mnum3, cpos3)
