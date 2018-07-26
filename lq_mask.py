import logging, os, sys
import numpy  as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from lq_utils import write_fastq, parse_fastx
from lq_exec  import LqExec
from time     import sleep

class LqMask:

    def __init__(self, reads, path_to_sdust, work_dir, *, suffix=None, n_proc=5, logger=None):
        if suffix:
            self.suffix = "_" + suffix
        else:
            self.suffix = ""
        if logger:
            self.logger = logger
        else:
            self.logger = logging.getLogger(__name__)
            self.logger.setLevel(logging.INFO)
            sh = logging.StreamHandler()
            formatter = logging.Formatter('%(module)s:%(asctime)s:%(lineno)d:%(levelname)s:%(message)s')
            sh.setFormatter(formatter)
            self.logger.addHandler(sh)
        if not os.path.isdir(work_dir):
            os.makedirs(work_dir, exist_ok=True)
        self.reads  = reads
        self.n_proc = n_proc
        self.psdust = path_to_sdust
        self.wdir   = work_dir
        self.outf   = os.path.join(work_dir, "longqc_sdust" + self.suffix + ".txt")
        self.tin  = []
        self.tout = []

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
        self.logger.info("sdust output file %s was made." % self.outf)
        
        for tf in self.tin + self.tout:
            if os.path.exists(tf):
                try:
                    os.remove(tf)
                    self.logger.info("tmp file %s was removed." % tf)
                except (OSError, e):
                    self.logger.error("%s - %s." % (e.filename, e.strerror))
            else:
                self.logger.warning("tmp file %s does not exist. skip removal of this file.")

    def run_async_sdust(self):
        procs = []
        n_seqs = len(self.reads)

        if not os.path.isdir(os.path.join(self.wdir, "analysis")):
            os.makedirs(os.path.join(self.wdir, "analysis"), exist_ok=True)

        for i in np.arange(0, self.n_proc):
            s = int(i * n_seqs/self.n_proc)
            e = int((i+1) * n_seqs/self.n_proc)
            fp = os.path.join(self.wdir, "analysis", "tmp_"+str(i)+".fastq")
            self.tin.append(fp)
            self.logger.debug("Seqs from %d to %d" % (s, e))
            write_fastq(fp, self.reads[s:e])
            p = LqExec(self.psdust, logger=self.logger)
            fpo = os.path.join(self.wdir, "analysis", "tmp_" + str(i) + self.suffix + ".txt")
            self.tout.append(fpo)
            p.exec(fp, out=fpo)
            self.logger.info("sdust process %s started." % p.get_pid() )
            procs.append(p)
        while True:
            for p in procs:
                if p.get_poll() is not None:
                    self.logger.info("sdust process %s terminated." % p.get_pid() )
                    procs.remove(p)
            self.logger.info("Calculating low complexity region...")

            if len(procs) == 0:
                break
            else:
                sleep(5)
        self.logger.info("Calculation finished.")
        self._concat_and_remove_tfiles()

    def get_outfile_path(self):
        return self.outf

# test
if __name__ == "__main__":

    # test
    lm = LqMask(reads_obj, bin_path, dir_path, n_proc=10)
    lm.run_async_sdust()
    sdust_outf = lm.get_outfile_path()
    print(sdust_outf)
    lm.plot_masked_fraction()
