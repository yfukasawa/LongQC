import sys, os, pysam, shutil
import numpy  as np

# https://stackoverflow.com/questions/5574702/how-to-print-to-stderr-in-python
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

#https://stackoverflow.com/questions/1868714/how-do-i-copy-an-entire-directory-of-files-into-an-existing-directory-using-pyth
def copytree(src, dst, symlinks=False, ignore=None):
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            shutil.copytree(s, d, symlinks, ignore)
        else:
            shutil.copy2(s, d)

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

def open_seq(fn):
    file_code = guess_format(fn)
    if file_code == 0:
        (reads, n_seqs, n_bases) =  parse_bam(fn, is_sequel=True)
        return (file_code, reads, n_seqs, n_bases)
    elif file_code == 1:
        eprint("Sorry. SAM is not supported.")
        return -1
    elif file_code == 2:
        (reads, n_seqs, n_bases) = parse_fastx(fn)
        return (file_code, reads, n_seqs, n_bases)
    elif file_code == 3:
        (reads, n_seqs, n_bases) = parse_fastx(fn)
        return (file_code, reads, n_seqs, n_bases)
    else:
        eprint("Sorry. the input file format is unknown and not supported.")
        return -1

# 0->bam, 1->sam, 2->fastq, 3->fasta, -1->error
def guess_format(fn):
    try:
        fh = open(fn, 'rb')
    except:
        eprint("Error: cannot open %s" % fn)

    try:
        majic = os.read(fh.fileno(), 4)
    except:
        eprint("Error: cannot read %s" % fn)

    # pybam and/or biopython way
    if majic == 'BAM\1':
        return 0
        fh.close()
    elif majic == b'\x1f\x8b\x08\x04':
        # compressed bam
        return 0
        fh.close()

    fh.close()

    try:
        fh = open(fn, 'r')
    except:
        eprint("Error: cannot open %s" % fn)

    # assume sam, fastx
    at_line_cnt = 0
    for line in fh:
        if line[0] == '@':
            at_line_cnt += 1
            continue
        elif at_line_cnt > 0:
            if at_line_cnt > 1:
                # header of sam
                fh.close()
                return 1
            cn = len(line.split("\t"))
            if cn == 11:
                fh.close()
                return 1
            at_line_cnt = 0
            # fastq
            fh.close()
            return 2
        elif line[0] == '>' and at_line_cnt == 0:
            # fasta
            fh.close()
            return 3
        else:
            cn = len(line.split("\t"))
            if cn == 11:
                fh.close()
                return 1
            at_line_cnt = 0
            continue

    # something else
    fh.close()
    return -1

# this is basically for pbbam, where there is no QV, and no support from minimap2.
# is_sequel option should be True, otherwise, this'll be very slow.
# https://github.com/lh3/minimap2/issues/159
def parse_bam(fn, *, is_sequel=True):
    reads   = []
    n_seqs  = 0
    n_bases = 0

    # basically for sequel
    input_bam = pysam.AlignmentFile(fn, 'rb', check_sq=False)
    for e in input_bam:
        n_seqs  += 1
        n_bases += len(e.query_sequence)
            
        if is_sequel:
            qual_33 = '!' * e.query_length
        else:
            qual_33 = "".join([chr(q+33) for q in e.query_qualities])
    
        reads.append( [e.query_name, e.query_sequence, qual_33] )

    return (reads, n_seqs, n_bases)

def parse_fastx(fn):
    reads   = []
    n_seqs  = 0
    n_bases = 0

    with pysam.FastxFile(fn) as f:
        for e in f:
            if e.quality:
                reads.append( [e.name, e.sequence, e.quality] )
            else:
                reads.append( [e.name, e.sequence] )
            n_seqs  += 1
            n_bases += len(e.sequence)

    return (reads, n_seqs, n_bases)

# deprecated.
def __parse_fastq(fn):
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
            reads.append( [name, seq, qual] )
    return (reads, n_seqs, n_bases)

def get_Qx_bases(reads, threshold=10):
    _t  = threshold + 33
    num = 0

    if len(reads[0]) < 3:
        # fasta case
        return num

    for read in reads:
        for i in range(len(read[2])):
            if ord(read[2][i]) >= _t:
                num += 1

    return num

# deprecated. no support for multi line, but this is not practical.
def __parse_fasta(fn):
    reads = []
    n_seqs  = 0
    n_bases = 0
    with open(fn, 'r') as fq:
        for line in fq:
            n_seqs += 1
            name = line.strip()[1:]
            seq  = next(fq).strip()
            n_bases += len(seq)
            reads.append( [name, seq] )
    return (reads, n_seqs, n_bases)

def write_fastq(fn, reads):
    if(os.path.isfile(fn)):
        eprint("Error: the file %s already exists." % fn)
        return 1

    with open(fn, 'w') as fq:
        for r in reads:
            fq.write("@%s\n%s\n+\n%s\n" % tuple(r))

# follow the logic flow of seqtk
# 2018/4/30 added exclude list. This is an ad-hoc way, and violates sampling schema. but let's see.
def sample_random_fastq_list(reads, param, *, s_seed=7, elist=None):
    frac    = 0.
    num     = 0
    n_seqs  = 0
    s_n_seqs  = 0
    s_n_bases = 0
    s_reads   = []

    if param > 1.:
        num = param
        s_reads = [0] * num
    else:
        frac = param

    np.random.seed(seed=s_seed)
    for read in reads:
        if(n_seqs%100000 == 0):
            h = np.random.uniform(size = 100000)
        name = read[0]
        seq  = read[1]
        qual = read[2]
        if elist and name in elist:
            #print("%s is skipped." % name)
            continue
        n_seqs  += 1
        if num:
            if(n_seqs - 1 < num):
                d = n_seqs-1
            else:
                d = int(h[(n_seqs-1)%100000]*n_seqs)
            if(d < num):
                if s_reads[d]:
                    s_n_seqs  -= 1
                    s_n_bases -= len(s_reads[d][1])
                s_reads[d] = [name, seq, qual]
                s_n_seqs  += 1
                s_n_bases += len(seq)
        elif( h[(n_seqs-1)%100000] < frac):
            s_reads.append([name, seq, qual])
            s_n_seqs  += 1
            s_n_bases += len(seq)
    return (s_reads, s_n_seqs, s_n_bases)

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
            #n_bases += len(seq)
            if elist and name in elist:
                #eprint("%s is skipped." % name)
                continue
            n_seqs  += 1
            if num:
                if(n_seqs - 1 < num):
                    d = n_seqs-1
                else:
                    d = int(h[(n_seqs-1)%100000]*n_seqs)
                if(d < num):
                    if reads[d]:
                        s_n_seqs  -= 1
                        s_n_bases -= len(reads[d][1])
                    reads[d] = [name, seq, qual]
                    s_n_seqs  += 1
                    s_n_bases += len(seq)
            elif( h[(n_seqs-1)%100000] < frac):
                reads.append([name, seq, qual])
                s_n_seqs  += 1
                s_n_bases += len(seq)
        return (reads, s_n_seqs, s_n_bases)

# test
if __name__ == "__main__":
    # some test code here
    code = guess_format("/home/fukasay/analyses/ont/minimap2_mapping/Hai1D_albacore213_1dsq_map.sam")
    print(code)

