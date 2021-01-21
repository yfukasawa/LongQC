import sys, os, pysam, shutil
import base64
import lq_nanopore
import numpy  as np
import gzip

from logging import getLogger
logger = getLogger(__name__)

def enc_b64_str(file_path):
    with open(file_path, "rb") as f:
        enc = base64.b64encode(f.read())
    return enc.decode('utf-8')

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
    a = np.array(vals)
    t = a.sum()/2
    cnt = 0
    a[::-1].sort() # this cannot be evaluated properly in for-loop
    for x in a:
        cnt += x
        if cnt >= t: return x

def get_NXX(vals, target=90):
    a = np.array(vals)
    t = a.sum() * target/100
    if target < 0:
        return vals[0]
    elif target > 100:
        return vals[-1]
    cnt = 0
    a[::-1].sort() # this cannot be evaluated properly in for-loop
    for x in a:
        cnt += x
        if cnt >= t: return x

def open_seq_chunk(fn, file_code, is_upper=False, chunk_size=500*1024**2): # default 500MB
    #file_code = guess_format(fn)
    if file_code == 0:
        yield from parse_bam_chunk(fn, chunk_size, is_sequel=True, is_upper=is_upper)
    elif file_code == 4:
        yield from parse_fast5_chunk(fn, chunk_size, is_upper=is_upper)
    elif file_code == 1:
        logger.error("SAM is not supported.")
        return -1
    elif file_code == 2 or file_code == 3:
        yield from parse_fastx_chunk(fn, chunk_size, is_upper=is_upper)
    else:
        logger.error("The input file format is unknown and not supported yet.")
        return -1

def open_seq(fn):
    file_code = guess_format(fn)
    if file_code == 0:
        (reads, n_seqs, n_bases) =  parse_bam(fn, is_sequel=True)
        return (file_code, reads, n_seqs, n_bases)
    elif file_code == 1:
        logger.error("Sorry. SAM is not supported.")
        return -1
    elif file_code == 2:
        (reads, n_seqs, n_bases) = parse_fastx(fn)
        return (file_code, reads, n_seqs, n_bases)
    elif file_code == 3:
        (reads, n_seqs, n_bases) = parse_fastx(fn)
        return (file_code, reads, n_seqs, n_bases)
    else:
        logger.error("Sorry. the input file format is unknown and not supported.")
        return -1

# 0->bam, 1->sam, 2->fastq(.gz), 3->fasta(.gz), 4->fast5 (multi), -1->error
def guess_format(fn):

    # assume fast5 is given in a dir. 
    if os.path.isdir(fn):
        logger.info("not a file but a direcory %s is given. looking for fast5 files.." % fn)
        for f in os.listdir(fn):
            if f.endswith(".fast5"):
                f5 = lq_nanopore.open_fast5(os.path.join(fn, f))
                if '/UniqueGlobalKey' in f5:
                    logger.error("single read fast5 is included? it's not supported for sampleqc.")
                    return -1
                return 4

        logger.error("no fast5 is found.")
        return -1

    try:
        fh = open(fn, 'rb')
    except:
        logger.error("cannot open %s" % fn)

    try:
        majic = os.read(fh.fileno(), 4)
    except:
        logger.error("cannot read %s" % fn)

    # pybam and/or biopython way
    if majic == 'BAM\1':
        fh.close()
        logger.debug("%s is an uncompressed BAM." % fn)
        return 0
    elif b'\x1f\x8b' in majic:
        # YF memo: 1f 8b 08 04 code can exist in fq.gz either.
        # changed the logic.
        fh.close()
        with gzip.open(fn, 'rb') as f:
            l = f.read(4)
            if "BAM" in l.decode(): # this should be 'BAM\x01'
                logger.debug("%s is a compressed BAM." % fn)
                return 0
            else:
                return __guess_sam_fastx(fn, isgzip=True)
    else:
        fh.close()

    return __guess_sam_fastx(fn, isgzip=False)

# guessing format is sam, fastq, fasta or something else.
def __guess_sam_fastx(fn, isgzip=False):
    if isgzip:
        try:
            fh = gzip.open(fn, 'rt')
        except:
            logger.error("cannot open a gzip file %s" % fn)
            sys.exit(1)
    else:
        try:
            fh = open(fn, 'r')
        except:
            logger.error("cannot open %s" % fn)
            sys.exit(1)

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


def parse_fast5_chunk(dn, cs, is_upper=False):
    reads   = []
    n_seqs  = 0
    n_bases = 0
    size = 0
    f5s = [os.path.join(dn, f) for f in os.listdir(dn) if f.endswith(".fast5")]
    for f5 in f5s:
        f5h = lq_nanopore.open_fast5(f5)
        top = lq_nanopore.list_toplevel(f5h)
        for k in top:
            if not k.startswith('read_'):
                continue
            fastq = lq_nanopore.get_fastq_from_multi_fast5(f5h, k).splitlines()
            name = fastq[0].split(" ")[0]
            if is_upper:
                reads.append( [name, fastq[1].upper(), fastq[3]] )
            else:
                reads.append( [name, fastq[1], fastq[3]] )
            size += sys.getsizeof(name) + sys.getsizeof(fastq[1]) + sys.getsizeof(fastq[3])
            n_bases += len(fastq[1])
            n_seqs += 1
            if size >= cs:
                yield (reads, n_seqs, n_bases)
                size  = 0
                reads = []
    yield (reads, n_seqs, n_bases)

def parse_bam_chunk(fn, cs, is_sequel=True, is_upper=False):
    reads   = []
    n_seqs  = 0
    n_bases = 0
    size = 0
    # basically for sequel
    input_bam = pysam.AlignmentFile(fn, 'rb', check_sq=False)
    for e in input_bam:
        n_seqs  += 1
        n_bases += len(e.query_sequence)
        if is_sequel:
            qual_33 = '!' * e.query_length
        else:
            qual_33 = "".join([chr(q+33) for q in e.query_qualities])
        if is_upper:
            reads.append( [e.query_name, e.query_sequence.upper(), qual_33] )
        else:
            reads.append( [e.query_name, e.query_sequence, qual_33] )
        size += sys.getsizeof(e.query_name) + sys.getsizeof(e.query_sequence) + sys.getsizeof(qual_33)
        if size >= cs:
            yield (reads, n_seqs, n_bases)
            size  = 0
            reads = []
    yield (reads, n_seqs, n_bases)

def parse_fastx_chunk(fn, cs, is_upper=False):
    reads   = []
    n_seqs  = 0
    n_bases = 0
    size = 0
    with pysam.FastxFile(fn) as f:
        for e in f:
            a = []
            if e.quality:
                if is_upper:
                    reads.append( [e.name, e.sequence.upper(), e.quality] )
                else:
                    reads.append( [e.name, e.sequence, e.quality] )
                size += sys.getsizeof(e.name) + sys.getsizeof(e.sequence) + sys.getsizeof(e.quality)
            else:
                if is_upper:
                    reads.append( [e.name, e.sequence.upper(), "!"*len(e.sequence)] )
                else:
                    reads.append( [e.name, e.sequence, "!"*len(e.sequence)] )
                size += sys.getsizeof(e.name) + sys.getsizeof(e.sequence) + sys.getsizeof("!"*len(e.sequence))
            n_seqs  += 1
            n_bases += len(e.sequence)
            if size >= cs:
                yield (reads, n_seqs, n_bases)
                size  = 0
                reads = []
    yield (reads, n_seqs, n_bases)

def parse_fastx(fn):
    reads   = []
    n_seqs  = 0
    n_bases = 0

    with pysam.FastxFile(fn) as f:
        for e in f:
            if e.quality:
                reads.append( [e.name, e.sequence, e.quality] )
            else:
                reads.append( [e.name, e.sequence, "!"*len(e.sequence)] )
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

def write_fastq(fn, reads, is_chunk=False):
    if(is_chunk==False and os.path.isfile(fn)):
        logger.error("the file %s already exists." % fn)
        return None

    try:
        reads[0]
    except IndexError:
        logger.error("No read to be output")
        return None

    # due to chunking, write mode was changed to a.
    mode = 'a' if is_chunk else 'w'
    with open(fn, mode) as fq:
        for r in reads:
            fq.write("@%s\n%s\n+\n%s\n" % tuple(r))

    return True

def subsample_from_chunk(chunk, cum_n_seq, s_reads, param, s_seed=7, elist=None):
    frac    = 0.
    num     = 0
    n_seqs  = cum_n_seq
    k = 0

    if param >= 1.:
        num = param
        if not s_reads:
            logger.info("list for subsample is not initialized. Initializing now.")
            s_reads = [0] * num
    else:
        frac = param
        a = []

    np.random.seed(seed=s_seed)
    h = np.random.uniform(size = len(chunk)+1)
        
    for read in chunk:
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
                d = int(h[k]*n_seqs)
            if(d < num):
                s_reads[d] = [name, seq, qual]
        elif( h[k] < frac):
            a.append([name, seq, qual])
        k += 1

    if num:
        return s_reads
    else:
        return s_reads + a

# follow the logic flow of seqtk
# 2018/4/30 added exclude list. This is an ad-hoc way, and violates sampling schema. but let's see.
def sample_random_fastq_list(reads, param, *, s_seed=7, elist=None):
    frac    = 0.
    num     = 0
    n_seqs  = 0
    s_n_seqs  = 0
    s_n_bases = 0
    s_reads   = []

    if param >= 1.:
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

    if num and n_seqs < num:
        s_reads = s_reads[:n_seqs]
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
                #logger.error("%s is skipped." % name)
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
    #code = guess_format("/home/fukasay/analyses/ont/minimap2_mapping/Hai1D_albacore213_1dsq_map.sam")
    #print(code)

    #fn = "/home/fukasay/rawdata/pb/rs2_ecoli_pacbio_official/ecoli_pacbio_p6c4.subreads.bam"
    fn = "/home/fukasay/python_codes/longQC/test_data"
    if guess_format(fn) == 4:
        for (reads, n_seqs, n_bases) in parse_fast5_chunk(fn, 0.5*1024**3):
            print(n_seqs, n_bases)
    #(c, reads, n_seqs, n_bases) = open_seq(fn)
    #print(len(reads), n_seqs, n_bases)
    #for (reads, n_seqs, n_bases) in open_seq_chunk(fn, guess_format(fn), 0.25*1024**3):
    #    print(len(reads), n_seqs, n_bases)
