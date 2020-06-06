import random
import array
import sys
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot  as plt
import matplotlib.mlab    as ml
import numpy              as np
from scipy.stats import gaussian_kde
from lq_utils    import guess_format, open_seq_chunk, open_seq

from logging import getLogger
logger = getLogger(__name__)

class LqGC:
    def __init__(self, chunk_size=150):
        self.chunk_size = chunk_size
        self.r_frac  = array.array('f')
        self.c_frac  = array.array('f')
        self.r_tot = 0
        self.c_tot = 0
        self.r_gc_tot = 0
        self.c_gc_tot = 0

    def calc_read_and_chunk_gc_frac(self, reads, samp_rate=0.2):
        # data should be capitalized
        for r in reads:
            s    = r[1]
            l    = len(s)
            self.r_tot += l
            gc_n = 0
            gc_n += s.count('G')
            gc_n += s.count('C')
            self.r_frac.append(gc_n/l)
            self.r_gc_tot += gc_n

            indices = np.random.choice(l, int(float(1/self.chunk_size)*l*samp_rate), replace=False)
            for i in indices:
                if i + self.chunk_size -1 > l:
                    break
                j = i + self.chunk_size
                cgc_n = 0
                cgc_n += s.count('G', i, j)
                cgc_n += s.count('C', i, j)
                self.c_frac.append(float(cgc_n)/self.chunk_size)
                self.c_gc_tot += cgc_n
                self.c_tot += self.chunk_size

    def plot_unmasked_gc_frac(self, fp=None, b_width = 0.02):

        dens_read = None

        logger.info("Mean GC composition: %.3f" % float(self.r_gc_tot/self.r_tot) )
        rtn_list = [np.mean(self.r_frac), np.std(self.r_frac)]
        plt.hist(self.r_frac, alpha=0.3, bins=np.arange(min(self.r_frac), max(self.r_frac) + b_width, b_width), color='blue', density=True)
        if len(self.r_frac) > 1:
            dens_read = gaussian_kde(self.r_frac)
        logger.info("Kernel density estimation done for read GC composition")

        plt.hist(self.c_frac, alpha=0.3, bins=np.arange(min(self.c_frac), max(self.c_frac) + b_width, b_width), color='red', density=True)
        logger.debug("Length of chunk array: %d" % len(self.c_frac))
        dens_chunk = gaussian_kde(self.c_frac)
        logger.info("Kernel density estimation done for chunked read GC composition")
        plt.grid(True)
        xs = np.linspace(0,1.0,50)
        if dens_read:
            plt.plot(xs, dens_read(xs), label="GC fraction read")
        plt.plot(xs, dens_chunk(xs), label="GC fraction of chunked read "+"("+str(self.chunk_size)+ "bp)")
        logger.debug("mean %f, stdev %f" % (np.mean(self.c_frac), np.std(self.c_frac)))
        #plt.plot(xs, ml.normpdf(xs, np.mean(res[0]), np.std(res[0])), label="GC fraction of chunked read "+"("+str(CHUNK_SIZE)+ "bp)")
        plt.xlabel('GC fraction')
        plt.ylabel('Probability density')

        plt.legend(bbox_to_anchor=(1,1), loc='upper right', borderaxespad=1)
        if fp:
            plt.savefig(fp, bbox_inches="tight", transparent=True)
        else:
            plt.show()
        plt.close()

        # bw param (https://stackoverflow.com/questions/23630515/getting-bandwidth-used-by-scipys-gaussian-kde-function)
        # to reproduce plots, we need not only bw but also data points. mmm, tentatively skip.
        #rtn_list.append(dens_read.covariance_factor()*np.array(self.r_frac).std())
        #rtn_list.append(dens_chunk.covariance_factor()*np.array(self.c_frac).std())
        return rtn_list

    """
    # deprecated
    def calc_masked_chunk_read_gc_frac(self, reads, chunk_size=300, th=0.2):
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

    def calc_chunk_read_gc_frac(self, reads, samp_rate=0.2, chunk_size=300):
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
        return (gc_f, tot_gc, tot)

    def calc_masked_read_gc_frac(self, reads):
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

    def calc_read_gc_frac(self, reads):
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

    def plot_unmasked_gc_frac_old(self, reads, fp=None, CHUNK_SIZE = 150, b_width = 0.02):

        dens_read = None

        res = self.calc_read_gc_frac(reads)
        logger.info("Mean GC composition: %.3f" % float(res[1]/res[2]) )
        rtn_list = [np.mean(res[0]), np.std(res[0])]
        plt.hist(res[0], alpha=0.3, bins=np.arange(min(res[0]), max(res[0]) + b_width, b_width), color='blue', normed=True)
        if len(res[0]) > 1:
            dens_read = gaussian_kde(res[0])
        logger.info("Kernel density estimation done for read GC composition")

        res = self.calc_chunk_read_gc_frac(reads, chunk_size=CHUNK_SIZE)
        plt.hist(res[0], alpha=0.3, bins=np.arange(min(res[0]), max(res[0]) + b_width, b_width), color='red', normed=True)
        logger.debug("Length of chunk array: %d" % len(res[0]))
        dens_chunk = gaussian_kde(res[0])
        logger.info("Kernel density estimation done for chunked read GC composition")
        plt.grid(True)
        xs = np.linspace(0,1.0,200)
        if dens_read:
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
    """

# test
if __name__ == "__main__":
    fn   = sys.argv[1]
    outf = sys.argv[2]

    # test
    lg = LqGC(chunk_size=150)
    file_format_code = guess_format(fn)
    chunk_n = 0
    for (reads, n_seqs, n_bases) in open_seq_chunk(fn, file_format_code, 0.5*1024**3):
        print("Computation of the GC fraction started for a chunk %d" % chunk_n)
        lg.calc_read_and_chunk_gc_frac(reads)
        chunk_n += 1
    plt.rcParams['figure.figsize'] = (7, 7)
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    print(lg.plot_unmasked_gc_frac(fp=outf))

