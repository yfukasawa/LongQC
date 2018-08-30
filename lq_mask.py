import os, sys
import multiprocessing as mp
import subprocess
import numpy  as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from lq_utils import write_fastq, open_seq_chunk, guess_format
from lq_exec  import LqExec
from time     import sleep

from logging import getLogger
logger = getLogger(__name__)

def _sdust(psdust, fin, fout):
    f = open(fout, "w")
    completed_process = subprocess.run([psdust, fin], check=True, stdout=f)
    if completed_process.stderr:
        return completed_process.stderr.decode('utf-8')
    else:
        return None

class LqMask:
    def __init__(self, path_to_sdust, work_dir, reads=None, suffix=None, max_n_proc=5):
        if suffix:
            self.suffix = "_" + suffix
        else:
            self.suffix = ""
        if not os.path.isdir(work_dir):
            os.makedirs(work_dir, exist_ok=True)
        if reads:
            self.reads  = reads
        self.n_proc = max_n_proc
        self.psdust = path_to_sdust
        self.wdir   = work_dir
        self.outf   = os.path.join(work_dir, "longqc_sdust" + self.suffix + ".txt")
        self.tin  = []
        self.tout = []
        self.pool = mp.Pool(self.n_proc)

    def plot_qscore_dist(self, df, column_qv, column_length, *, fp=None, platform='ont', interval=3000):
        if platform == 'ont':
            mid_threshold = 7 # ont
        else:
            mid_threshold = 8 # pb
        df['Binned read length'] = np.floor(df[column_length].values/interval)
        df.boxplot(column=column_qv, by='Binned read length', sym='+', rot=90, figsize=(2*int(max(df['Binned read length'])/5+0.5), 4.8))
        plt.grid(True)
        xmin, xmax = plt.gca().get_xlim()
        ymin, ymax = plt.gca().get_ylim()
        plt.xticks(np.arange(xmax+1), [int(i) for i in np.arange(xmax+1)*interval])
        plt.axhspan(0,  mid_threshold, facecolor='red', alpha=0.1)
        #plt.axhspan(5,  mid_threshold, facecolor='yellow', alpha=0.1)
        plt.axhspan(mid_threshold, ymax, facecolor='green', alpha=0.1)
        #plt.boxplot(df[5].values[np.where(df[4] == 0.0)])
        plt.ylim(0, ymax)
        plt.ylabel('Averaged QV')
        plt.title("")
        plt.suptitle("")
        if fp:
            plt.savefig(fp, bbox_inches="tight")
        else:
            plt.show()
        plt.close()

    def plot_masked_fraction(self, fp=None):
        self.df = pd.read_table(self.outf, sep='\t', header=None)
        plt.grid(True)
        plt.hist(self.df[3], alpha=0.2, bins=np.arange(0, 1.0, 0.01), color='red')
        plt.xlim(0, 1.0)
        plt.xlabel('Low complexity fraction')
        plt.ylabel('Frequency')
        if fp:
            plt.savefig(fp, bbox_inches="tight")
        else:
            plt.show()
        plt.close()

    def _concat_and_remove_tfiles(self):
        with open(self.outf, 'w') as out:
            for tf in self.tout:
                with open(tf, 'r') as t:
                    for l in t:
                        out.write(l)
        logger.info("sdust output file %s was made." % self.outf)
        
        for tf in self.tin + self.tout:
            if os.path.exists(tf):
                try:
                    os.remove(tf)
                    logger.info("tmp file %s was removed." % tf)
                except (OSError, e):
                    logger.error("%s - %s." % (e.filename, e.strerror))
            else:
                logger.warning("tmp file %s does not exist. skip removal of this file.")

    # for multiple call case like chunking
    def submit_sdust(self, reads, chunk_n):
        if not os.path.isdir(os.path.join(self.wdir, "analysis")):
            logger.info("A new dir was made: %s" % os.path.join(self.wdir, "analysis"))
            os.makedirs(os.path.join(self.wdir, "analysis"), exist_ok=True)

        fpi = os.path.join(self.wdir, "analysis", "tmp_" + str(chunk_n) + ".fastq")
        self.tin.append(fpi)
        fpo = os.path.join(self.wdir, "analysis", "tmp_" + str(chunk_n) + self.suffix + ".txt")
        self.tout.append(fpo)
        write_fastq(fpi, reads)
        self.pool.apply_async(_sdust, args=(self.psdust, fpi, fpo))
        logger.info("New job was submitted: in->%s, out->%s" % (fpi, fpo))

    def close_pool(self):
        logger.info("Waiting completion of all of jobs...")
        self.pool.close()
        self.pool.join()
        logger.info("sdust jobs finished.")
        self._concat_and_remove_tfiles()

    # for a single call case
    def run_async_sdust(self):
        procs = []
        if self.reads:
            n_seqs = len(self.reads)
        else:
            logger.error("No read is given for analysis.")
            sys.exit(1)

        if not os.path.isdir(os.path.join(self.wdir, "analysis")):
            os.makedirs(os.path.join(self.wdir, "analysis"), exist_ok=True)

        for i in np.arange(0, self.n_proc):
            s = int(i * n_seqs/self.n_proc)
            e = int((i+1) * n_seqs/self.n_proc)
            fp = os.path.join(self.wdir, "analysis", "tmp_"+str(i)+".fastq")
            self.tin.append(fp)
            logger.debug("Seqs from %d to %d" % (s, e))
            write_fastq(fp, self.reads[s:e])
            p = LqExec(self.psdust)
            fpo = os.path.join(self.wdir, "analysis", "tmp_" + str(i) + self.suffix + ".txt")
            self.tout.append(fpo)
            p.exec(fp, out=fpo)
            logger.info("sdust process %s started." % p.get_pid() )
            procs.append(p)
        while True:
            for p in procs:
                if p.get_poll() is not None:
                    logger.info("sdust process %s terminated." % p.get_pid() )
                    procs.remove(p)
            logger.info("Calculating low complexity region...")

            if len(procs) == 0:
                break
            else:
                sleep(5)
        logger.info("Calculation finished.")
        self._concat_and_remove_tfiles()

    def get_outfile_path(self):
        return self.outf

# test
if __name__ == "__main__":
    # test
    
    lm = LqMask("sdust", "./")
    chunk_n = 0
    fn = sys.argv[1]
    file_code = guess_format(fn)
    for (reads, n_seqs, n_bases) in open_seq_chunk(fn, file_code, chunk_size=float(sys.argv[2])*1024**3, is_upper=True):
        lm.submit_sdust(reads, chunk_n)
        chunk_n += 1
    lm.close_pool()
    lm.plot_masked_fraction("./masked_frac.png")
