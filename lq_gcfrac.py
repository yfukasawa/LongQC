import random
import logging
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot  as plt
import matplotlib.mlab    as ml
import numpy              as np

from scipy.stats import gaussian_kde
from lq_utils    import parse_fastx

def calc_masked_chunk_read_gc_frac(reads, chunk_size=300, th=0.2):
    tot     = 0
    tot_gc  = 0

    gc_f    = []

    for r in reads:
        s = r[1]
        l = len(s)
        tot += l

        chunk = [s[i:i+chunk_size] for i in range(0, l, chunk_size)]
        for c in chunk:
            gc_n = 0
            gc_n += c.count('G')
            gc_n += c.count('C')
            tot_gc += gc_n
            if c.count('N')/chunk_size > th:
                continue
            gc_f.append(gc_n/chunk_size)

    return (gc_f, tot_gc, tot)


def calc_chunk_read_gc_frac(reads, samp_rate=0.2, chunk_size=300):
    tot     = 0
    tot_gc  = 0

    gc_f    = []
    
    #original codes generates too many sampling point, so randomly sample the points.
    for r in reads:
        s    = r[1]
        l    = len(s)

        indices = np.random.choice(len(s), int(float(1/chunk_size)*len(s)*samp_rate), replace=False )
        for i in indices:
            if i + chunk_size -1 > len(s):
                break
            j = i+ chunk_size
            gc_n = 0
            gc_n += s.count('G', i, j)
            gc_n += s.count('C', i, j)
            gc_n += s.count('g', i, j)
            gc_n += s.count('c', i, j)
            gc_f.append(float(gc_n)/float(chunk_size))
            tot_gc += gc_n
            tot += chunk_size

        """
        #chunk  = [s[i:i+chunk_size] for i in range(0, len(s), chunk_size)]
        for c in chunk:
            if len(c) < chunk_size:
                continue
            gc_n = 0
            gc_n += c.count('G')
            gc_n += c.count('C')
            gc_n += c.count('g')
            gc_n += c.count('c')
            gc_f.append(float(gc_n)/float(chunk_size))
            tot_gc += gc_n
            tot += chunk_size
        """

    return (gc_f, tot_gc, tot)

def calc_masked_read_gc_frac(reads):
    tot     = 0
    tot_gc  = 0

    gc_f    = []
    #reads can be fastq or fasta. second item of each tuple must be seq.
    for r in reads:
        s    = r[1]
        l    = len(s)
        _l   = l - s.count('N')
        tot += l
        
        gc_n = 0
        gc_n += s.count('G')
        gc_n += s.count('C')
        tot_gc += gc_n
        if s.count('N')/l > 0.9:
            continue
        gc_f.append(gc_n/_l)


    return (gc_f, tot_gc, tot)


def calc_read_gc_frac(reads):
    tot     = 0
    tot_gc  = 0

    gc_f    = []
    #reads can be fastq or fasta. second item of each tuple must be seq.
    for r in reads:
        s    = r[1]
        l    = len(s)
        tot += l
        
        gc_n = 0
        gc_n += s.count('G')
        gc_n += s.count('C')
        gc_n += s.count('g')
        gc_n += s.count('c')

        # below code might add unexpected bias. let's keep it simple.
        #amb_n  = 0
        #amb_n += s.count('N') + s.count('n')
        #gc_n  += amb_n * 0.25 

        gc_f.append(gc_n/l)
        tot_gc += gc_n


    return (gc_f, tot_gc, tot)


def plot_unmasked_gc_frac(reads, *, fp=None, logger=None, CHUNK_SIZE = 150, b_width = 0.02):

    if not logger:
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.DEBUG)
        sh = logging.StreamHandler()
        formatter = logging.Formatter('%(module)s:%(asctime)s:%(lineno)d:%(levelname)s:%(message)s')
        sh.setFormatter(formatter)
        logger.addHandler(sh)

    res = calc_read_gc_frac(reads)
    logger.info("Mean GC composition: %.3f" % float(res[1]/res[2]) )
    rtn_list = [np.mean(res[0]), np.std(res[0])]
    plt.hist(res[0], alpha=0.3, bins=np.arange(min(res[0]), max(res[0]) + b_width, b_width), color='blue', normed=True)
    dens_read = gaussian_kde(res[0])
    logger.info("Kernel density estimation done for read GC composition")

    res = calc_chunk_read_gc_frac(reads, chunk_size=CHUNK_SIZE)

    plt.hist(res[0], alpha=0.3, bins=np.arange(min(res[0]), max(res[0]) + b_width, b_width), color='red', normed=True)
    dens_chunk = gaussian_kde(res[0])
    logger.info("Kernel density estimation done for chunked read GC composition")
    plt.grid(True)
    xs = np.linspace(0,1.0,200)
    plt.plot(xs, dens_read(xs), label="GC fraction read")
    plt.plot(xs, dens_chunk(xs), label="GC fraction of chunked read "+"("+str(CHUNK_SIZE)+ "bp)")
    logger.debug("mean %f, stdev %f" % (np.mean(res[0]), np.std(res[0])))
    #plt.plot(xs, ml.normpdf(xs, np.mean(res[0]), np.std(res[0])), label="GC fraction of chunked read "+"("+str(CHUNK_SIZE)+ "bp)")
    plt.xlabel('GC fraction')
    plt.ylabel('Probability density')

    plt.legend(bbox_to_anchor=(1,1), loc='upper right', borderaxespad=1)
    if fp:
        plt.savefig(fp, bbox_inches="tight")
    else:
        plt.show()
    plt.close()

    return rtn_list


# test
if __name__ == "__main__":

    is_masked  = False

    reads, n_seqs, n_bases = parse_fastx(path_to_reads)
    print('reads were loaded.\n')

    plot_unmasked_gc_frac(reads)

    if is_masked:
        res = calc_masked_read_gc_frac(reads)
        print(res[1], res[2], "-->", res[1]/res[2])   

    if is_masked:
        res = calc_masked_chunk_read_gc_frac(reads, chunk_size=CHUNK_SIZE)
        print(res[1], res[2], "-->", res[1]/res[2])

