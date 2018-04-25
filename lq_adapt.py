import os, logging
import bindings.python.edlib

def _cut(reads, adp, th, r, lengths=[]):
    iden_max = -1

    for read in reads:
        lengths.append( len(read[1]) )
        #result5 = bindings.python.edlib.align(, read[1][:150], mode="HW", task='path')
        result = bindings.python.edlib.align(adp, read[1][:r], mode="HW", task='path')
        if result['identity'] >th:
            #print(read[0], result['identity'], result['cigar'])
            if result['identity'] > iden_max:
                iden_max = result['identity']

    return iden_max


def cut_adapter(reads, lengths, logger, adp_t, adp_b = None, th=0.75, length=150):
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
        iden_max5 = _cut(reads, adp_t, th, length, lengths)

    if adp_b:
        iden_max3 = _cut(reads, adp_t, th, length)
