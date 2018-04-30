import re
import sys
import os

import numpy  as np

def rgb(r,g,b):

    return [ r/255, g/255, b/255 ]

def get_N50(vals):
    t = np.array(vals).sum()/2
    cnt = 0
    for x in np.sort(vals):
        cnt += x
        if cnt >= t: return x

def get_NXX(vals, target=90):
    t = np.array(vals).sum() * target/100
    if target < 0:
        return vals[0]
    elif target > 100:
        return vals[-1]
    cnt = 0
    for x in np.sort(vals):
        cnt += x
        if cnt >= t: return x

# https://stackoverflow.com/questions/5574702/how-to-print-to-stderr-in-python
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def parse_fastq(fn):
    reads   = []
    n_seqs  = 0
    n_bases = 0
    with open(fn, 'r') as fq:
        for line in fq:
            n_seqs += 1
            name = line.strip()[1:]
            seq  = next(fq).strip()
            n_bases += len(seq)
            next(fq)
            qual = next(fq).strip()
            reads.append( (name, seq, qual) )
    return (reads, n_seqs, n_bases)

def parse_fasta(fn):
    reads = []
    n_seqs  = 0
    n_bases = 0
    with open(fn, 'r') as fq:
        for line in fq:
            n_seqs += 1
            name = line.strip()[1:]
            seq  = next(fq).strip()
            n_bases += len(seq)
            reads.append( (name, seq) )
    return (reads, n_seqs, n_bases)

def write_fastq(fn, reads):
    if(os.path.isfile(fn)):
        eprint("Error: the file %s already exists." % fn)
        return 1

    with open(fn, 'w') as fq:
        for r in reads:
            fq.write("@%s\n%s\n+\n%s\n" % r)

# follow the logic flow of seqtk
# 2018/4/30 added exclude list. This is an ad-hoc way, and violates sampling schema. but let's see.
def sample_random_fastq(fn, param, *, s_seed=7, elist=None):
    frac    = 0.
    num     = 0
    n_seqs  = 0
    s_n_seqs  = 0
    s_n_bases = 0
    reads   = []

    if param > 1.:
        num = param
        reads = [0] * num
    else:
        frac = param

    np.random.seed(seed=s_seed)
    with open(fn, 'r') as fq:
        for line in fq:
            if(n_seqs%100000 == 0):
                h = np.random.uniform(size = 100000)
            name = line.strip()[1:].split(" ")[0]
            seq  = next(fq).strip()
            next(fq)
            qual = next(fq).strip()
            n_seqs  += 1
            #n_bases += len(seq)
            if name in elist:
                print("%s is skipped." % name)
                continue
            if num:
                if(n_seqs - 1 < num):
                    d = n_seqs-1
                else:
                    d = int(h[(n_seqs-1)%100000]*n_seqs)
                if(d < num):
                    if reads[d]:
                        s_n_seqs  -= 1
                        s_n_bases -= len(reads[d][1])
                    reads[d] = (name, seq, qual)
                    s_n_seqs  += 1
                    s_n_bases += len(seq)
            elif( h[(n_seqs-1)%100000] < frac):
                reads.append((name, seq, qual))
                s_n_seqs  += 1
                s_n_bases += len(seq)
        return (reads, s_n_seqs, s_n_bases)


# test
if __name__ == "__main__":
    #fq_tuple = parse_fastq("/home/fukasay/basecalled/ont/20171024_0704_20171024_ecoli_1D_test_albacore213/release.fastq")
    #print(fq_tuple[1:3])
    #write_fastq("/tmp/test2.fastq", fq_tuple[0])

    (reads, n_seqs, n_bases) = sample_random_fastq("/home/fukasay/rawdata/pb/r54224_20171114_141933/1_A01/m54224_171114_143939.subreads.fastq", 10000)
    print(n_seqs, n_bases)
    #write_fastq("/home/fukasay/basecalled/ont/20171218_1536_20171218_ECOLI_1DSQUARE_TEST/1dsq_sub10k.fastq", reads)
