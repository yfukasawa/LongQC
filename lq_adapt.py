import os, logging
import bindings.python.edlib

def _cut(reads, adp, th, r, lengths=[]):
    iden_max  = -1
    match_num = 0

    for read in reads:
        lengths.append( len(read[1]) )
        #result5 = bindings.python.edlib.align(, read[1][:150], mode="HW", task='path')
        result = bindings.python.edlib.align(adp, read[1][:r], mode="HW", task='path')
        if result['identity'] >th:
            print(read[0], result['identity'], result['cigar'])
            print(result)
            match_num += 1        
            if result['identity'] > iden_max:
                iden_max = result['identity']

    return (iden_max, match_num)

def cut_adapter(reads, lengths, adp_t, adp_b = None, logger=None, th=0.75, length=150):
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

    tot       = 0
    trim5     = 0
    trim3     = 0
    trimboth  = 0
    iden_max5 = 0.0
    iden_max3 = 0.0
    max3_seq  = ''

    if not adp_t and not adp_b:
        logger.WARNING()
        return None
    
    if adp_t:
        (iden_max5, mnum) = _cut(reads, adp_t, th, length, lengths)
        logger.info("max identity:%f The number of matched read: %d" % (iden_max5, mnum))

    if adp_b:
        (iden_max3, mnum) = _cut(reads, adp_t, th, length)
        logger.info("max identity:%f The number of matched read: %d" % (iden_max3, mnum))

    """
    print('traverse was finished. %d reads were processed.' % tot)
    print('Reads having adapter-like seq in both ends: %d / %d' % (trimboth, tot))
    print('Reads having adapter-like seq in 5\' ends: %d / %d' % (trim5 - trimboth, tot))
    print('Reads having adapter-like seq in 3\' ends: %d / %d' % (trim3 - trimboth, tot))
    print('Max identity in the reads: %.3f' % iden_max5)
    print('Max identity in the reads: %.3f' % iden_max3)
    print(max3_seq)
    """
