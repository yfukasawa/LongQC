import os, logging
import bindings.python.edlib

def _cutr(reads, adp, th, r, *, len_list=[]):
    iden_max  = -1
    match_num = 0
    cut_pos   = []
    if len(reads[0]) > 2:
        has_qual = True
    else:
        has_qual = False
    for read in reads:
        if not len_list:
            pass
        else:
            len_list.append( len(read[1]) )
        result = bindings.python.edlib.align(adp, read[1][-r:], mode="HW", task='path')
        if result['identity'] >th:
            #print(read[0], result['identity'], result['locations'], result['cigar'])
            start = result['locations'][0][0]
            cut_pos.append(start)
            match_num += 1        
            if result['identity'] > iden_max:
                iden_max = result['identity']
            read[1] = read[1][:start]
            if has_qual:
                read[2] = read[2][:start]
    return (iden_max, match_num, cut_pos)

def _cutf(reads, adp, th, r, *, len_list=[]):
    iden_max  = -1
    match_num = 0
    cut_pos   = []
    if len(reads[0]) > 2:
        has_qual = True
    else:
        has_qual = False
    for read in reads:
        if not len_list:
            pass
        else:
            len_list.append( len(read[1]) )
        result = bindings.python.edlib.align(adp, read[1][:r], mode="HW", task='path')
        if result['identity'] >th:
            #print(read[0], result['identity'], result['locations'], result['cigar'])
            end = result['locations'][0][1]
            cut_pos.append(end)
            match_num += 1        
            if result['identity'] > iden_max:
                iden_max = result['identity']
            read[1] = read[1][end+1:]
            if has_qual:
                read[2] = read[2][end+1:]
    return (iden_max, match_num, cut_pos)

def cut_adapter(reads, *, len_list=None, adp_t=None, adp_b = None, logger=None, th=0.75, length=150):
    """
    > sequel_adapter
    ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT
    > rs2_adapter
    ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT

    > 1d Y-t
    AATGTACTTCGTTCAGTTACGTATTGCT
    > 1d Y-b
    GCAATACGTAACTGAACGAAGT

    > rapid kit
    GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA
    """

    iden_max5 = 0.0
    iden_max3 = 0.0

    if not adp_t and not adp_b:
        logger.WARNING()
        return None

    if adp_t and adp_b:
        (iden_max5, mnum5, cpos5) = _cutf(reads, adp_t, th, length, len_list=len_list)
        (iden_max3, mnum3, cpos3) = _cutr(reads, adp_t, th, length)
        logger.info("Adapter Sequence: %s, max identity:%f and the number of trimmed reads: %d" % (adp_t, iden_max5, mnum5))
        logger.info("Adapter Sequence: %s, max identity:%f and the number of trimmed reads: %d" % (adp_b, iden_max3, mnum3))
        return ((iden_max5, mnum5, cpos5), (iden_max3, mnum3, cpos3))

    if adp_t and not adp_b:
        (iden_max5, mnum5, cpos5) = _cutf(reads, adp_t, th, length, len_list=len_list)
        logger.info("Adapter Sequence: %s, max identity:%f and the number of trimmed reads: %d" % (adp_t, iden_max5, mnum5))
        return (iden_max5, mnum, cpos)

    if not adp_t and adp_b:
        (iden_max3, mnum3, cpos3) = _cutr(reads, adp_t, th, length)
        logger.info("Adapter Sequence: %s, max identity:%f and the number of trimmed reads: %d" % (adp_b, iden_max3, mnum3))
        return (iden_max3, mnum3, cpos3)

    """
    print('traverse was finished. %d reads were processed.' % tot)
    print('Reads having adapter-like seq in both ends: %d / %d' % (trimboth, tot))
    print('Reads having adapter-like seq in 5\' ends: %d / %d' % (trim5 - trimboth, tot))
    print('Reads having adapter-like seq in 3\' ends: %d / %d' % (trim3 - trimboth, tot))
    print('Max identity in the reads: %.3f' % iden_max5)
    print('Max identity in the reads: %.3f' % iden_max3)
    print(max3_seq)
    """
