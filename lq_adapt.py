import os, re
import numpy as np
import edlib
from lq_utils    import open_seq

import logging
#from logging import getLogger
logger = logging.getLogger(__name__)

def _cutr(reads, adp, th, r, len_list=[]):
    iden_max  = -1
    match_num = 0
    cut_pos   = []
    skip_num  = 0
    repat = re.compile(r'(\d+)[DHIMNPSX=]{1}')
    if len(reads[0]) > 2:
        has_qual = True
    else:
        has_qual = False
    for read in reads:
        read_length = len(read[1])
        if not len_list:
            pass
        else:
            len_list.append( read_length )
        if read_length < 2*r:
            #logger.debug("%s is too short: %d bp. Skip." % (read[0], read_length))
            skip_num += 1
            continue
        result = edlib.align(adp, read[1][-r:], mode="HW", task='path')
        identity = 1.0 - float(result['editDistance']/np.sum([int(i) for i in repat.findall(result['cigar'])]))
        if identity >th:
            start = len(read[1]) - r + result['locations'][0][0]
            cut_pos.append(r-result['locations'][0][0])
            match_num += 1        
            if identity > iden_max:
                iden_max = identity
            read[1] = read[1][:start]
            if has_qual:
                read[2] = read[2][:start]

    logger.info("%d reads were skipped due to their short lengths." % skip_num)
    return (iden_max, match_num, cut_pos)

def _cutf(reads, adp, th, r, len_list=[]):
    iden_max  = -1
    match_num = 0
    cut_pos   = []
    skip_num  = 0
    repat = re.compile(r'(\d+)[DHIMNPSX=]{1}')
    if len(reads[0]) > 2:
        has_qual = True
    else:
        has_qual = False
    for read in reads:
        read_length = len(read[1])
        if not len_list:
            pass
        else:
            len_list.append( read_length )
        if read_length < 2*r:
            skip_num += 1
            #logger.debug("%s is too short: %d bp. Skip." % (read[0], read_length))
            continue
        result = edlib.align(adp, read[1][:r], mode="HW", task='path')
        identity = 1.0 - float(result['editDistance']/np.sum([int(i) for i in repat.findall(result['cigar'])]))
        if identity >th:
            end = result['locations'][0][1]
            cut_pos.append(end)
            match_num += 1        
            if identity > iden_max:
                iden_max = identity
            read[1] = read[1][end+1:]
            if has_qual:
                read[2] = read[2][end+1:]

    logger.info("%d reads were skipped due to their short lengths." % skip_num)
    return (iden_max, match_num, cut_pos)

def cut_adapter(reads, len_list=None, adp_t=None, adp_b = None, th=0.75, length=150):
    iden_max5 = 0.0
    iden_max3 = 0.0

    if not adp_t and not adp_b:
        logger.error("No adapter sequence is given.")
        return None

    if adp_t and adp_b:
        (iden_max5, mnum5, cpos5) = _cutf(reads, adp_t, th, length, len_list=len_list)
        logger.info("Adapter Sequence: %s, max identity:%f and the number of trimmed reads: %d" % (adp_t, iden_max5, mnum5))
        (iden_max3, mnum3, cpos3) = _cutr(reads, adp_b, th, length)
        logger.info("Adapter Sequence: %s, max identity:%f and the number of trimmed reads: %d" % (adp_b, iden_max3, mnum3))
        return ((iden_max5, mnum5, cpos5), (iden_max3, mnum3, cpos3))

    if adp_t and not adp_b:
        (iden_max5, mnum5, cpos5) = _cutf(reads, adp_t, th, length, len_list=len_list)
        logger.info("Adapter Sequence: %s, max identity:%f and the number of trimmed reads: %d" % (adp_t, iden_max5, mnum5))
        return (iden_max5, mnum5, cpos5)

    if not adp_t and adp_b:
        (iden_max3, mnum3, cpos3) = _cutr(reads, adp_b, th, length, len_list=len_list)
        logger.info("Adapter Sequence: %s, max identity:%f and the number of trimmed reads: %d" % (adp_b, iden_max3, mnum3))
        return (iden_max3, mnum3, cpos3)

# test
if __name__ == "__main__":
    seqf = "/path/to/seq"

    adp5 = "ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT"
    adp3 = "ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT"

    #logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    sh = logging.StreamHandler()
    formatter = logging.Formatter('%(module)s:%(asctime)s:%(lineno)d:%(levelname)s:%(message)s')
    sh.setFormatter(formatter)
    logger.addHandler(sh)

    logger.info("Input file: %s" % seqf)
    (file_format_code, reads, n_seqs, n_bases) = open_seq(seqf)
    logger.info("File was loaded.")
    (tuple_5, tuple_3) = cut_adapter(reads, adp_t=adp5, adp_b=adp3)
    
