import math, logging
import numpy             as np
import pandas            as pd
import scipy.stats       as st
import matplotlib.pyplot as plt

#from scipy.stats   import gumbel_r
from sklearn  import mixture
from operator import itemgetter

"""
#### deprecated
# df:  DataFrame of pandas
# k :  the parameter for the number of component in GMM.
def est_coverage_dist(df, column_i, k=2):
    m_f   = mixture.GaussianMixture(n_components=k).fit(df[column_i].values.reshape(-1,1),1)
    order = m_f.weights_/[e[0] for inner in m_f.covariances_ for e in inner] # w/\sigma
    m_i   = 1 if order[0] < order[1] else 0

    return (m_f.means_[m_i][0], m_f.covariances_[m_i][0][0])


#### deprecated
# df:  DataFrame of pandas
# k :  the parameter for the number of component in GMM.
# nbs: the number of bootstrap
# Empirically, we determined k=2, main component and background component.

# This returns relatively narrow region, and actually not so reliable...
def est_confidence_interval_mu(df, column_i, k=2, nbs=1000):

    #_means = np.zeros((k, nbs))
    #_covs  = np.zeros((k, nbs))

    _mus     = np.zeros(nbs)
    #_sigmas  = np.zeros(nbs)

    #df = df[df[column_i] != 0.0]

    # full model
    m_f   = mixture.GaussianMixture(n_components=k).fit(df[column_i].values.reshape(-1,1),1)
    order = m_f.weights_/[e[0] for inner in m_f.covariances_ for e in inner] # w/\sigma
    m_i   = 1 if order[0] < order[1] else 0
    #b_i   = m_i ^ 1

    #for i in range(0, nbs):
    i = 0
    while(i < nbs):
        x  = df[column_i].values[np.random.randint(df.shape[0], size=df.shape[0])]
        _m = mixture.GaussianMixture(n_components=k).fit(x.reshape(-1,1),1)
        if _m.converged_:
            _o   = _m.weights_/[e[0] for inner in _m.covariances_ for e in inner] # w/\sigma
            _m_i = 1 if _o[0] < _o[1] else 0
            _mus[i]    = _m.means_[_m_i][0]
            #_sigmas[i] = np.sqrt(_m.covariances_[_m_i][0][0])
            
            i += 1

    return st.t.interval(0.99, len(_mus)-1, loc=np.mean(_mus), scale=st.sem(_mus))
"""
class LqCoverage:

    UNMAPPED_FRACTION_THRESHOLD = 0.4
    UNMAPPED_FRACTION_PARAM_MIN = 0.05
    UNMAPPED_FRACTION_PARAM_MAX = 0.2
    QLENGTH_COLUMN  = 1
    N_MBASE_COLUMN  = 2
    COVERAGE_COLUMN = 4
    QV_COLUMN = 5

    def __init__(self, table_path, logger=None):
        self.df = pd.read_table(table_path, sep='\t', header=None)
        self.min_lambda = None
        self.max_lambda = None
        self.unmapped_frac_trimmed   = 0.0
        self.unmapped_frac_untrimmed = 0.0
        self.model = None
        self.model_main_comp = None
        self.mean_main = None
        self.cov_main  = None
        if logger:
            self.logger = logger
        else:
            self.logger = logging.getLogger(__name__)
            self.logger.setLevel(logging.INFO)
            sh = logging.StreamHandler()
            formatter = logging.Formatter('%(module)s:%(asctime)s:%(lineno)d:%(levelname)s:%(message)s')
            sh.setFormatter(formatter)
            self.logger.addHandler(sh)

        self.logger.info("Estimating coverage distribution..")
        self.__est_coverage()
        self.logger.info("Estimation of coverage distribution finished.")

    def __est_coverage(self):
        self.unmapped_frac_trimmed   = len(self.df[LqCoverage.COVERAGE_COLUMN].values[np.where(self.df[LqCoverage.COVERAGE_COLUMN] == 0.0)]) \
                                       / len(self.df[LqCoverage.QLENGTH_COLUMN])
        self.unmapped_frac_untrimmed = len(self.df[LqCoverage.N_MBASE_COLUMN].values[np.where(self.df[LqCoverage.N_MBASE_COLUMN] == 0.0)]) \
                                       / len(self.df[LqCoverage.QLENGTH_COLUMN])
        self.logger.info("Unmapped fraction: %.3f (naive), %.3f (coverage considered)" % (self.unmapped_frac_untrimmed, self.unmapped_frac_trimmed))

        if self.unmapped_frac_trimmed >= LqCoverage.UNMAPPED_FRACTION_THRESHOLD:
            self.logger.warning("The fraction of zero coverage read is high %.3f" % self.unmapped_frac_trimmed)
            self.min_lambda = -1 * math.log(self.unmapped_frac_trimmed - LqCoverage.UNMAPPED_FRACTION_PARAM_MIN)
            self.max_lambda = -1 * math.log(self.unmapped_frac_trimmed - LqCoverage.UNMAPPED_FRACTION_PARAM_MAX)
            range_str = str(self.min_lambda) + "-" + str(self.max_lambda)
            self.logger.warning("If and only if the data is healthy, very rough estimated coverage range is %s." % range_str)

        self.model_main_comp = self.__est_coverage_dist_gmm(k_i=2)
        self.model = self.model_main_comp[0]
        self.mean_main = self.model_main_comp[1]
        self.cov_main  = self.model_main_comp[2]

    def plot_coverage_dist(self, fp=None):
        gmm_x = np.linspace(0, self.mean_main+4*np.sqrt(self.cov_main), 5000)
        gmm_y = np.exp(self.model.score_samples(gmm_x.reshape(-1,1)))
        plt.grid(True)

        if self.min_lambda and self.max_lambda:
            pois_min = st.poisson(self.min_lambda)
            pois_max = st.poisson(self.max_lambda)
            pois_x = np.arange(int(self.mean_main+4 * np.sqrt(self.cov_main)) + 1)
            pois_y_min = pois_min.pmf(pois_x)
            pois_y_max = pois_max.pmf(pois_x)
            plt.xlim(0, 50)
            plt.plot(pois_x, pois_y_min, label="Fitted Model by Poisson model (" + "%.3f" % self.min_lambda + ")")
            plt.plot(pois_x, pois_y_max, label="Fitted Model by Poisson model (" + "%.3f" % self.max_lambda + ")")
        else:
            plt.plot(gmm_x, gmm_y, label="Fitted by Gaussian mixture model")

        if self.unmapped_frac_trimmed < LqCoverage.UNMAPPED_FRACTION_THRESHOLD:
            plt.hist(self.df[LqCoverage.COVERAGE_COLUMN],
                     alpha=0.2,
                     bins=np.arange(0,
                                    self.mean_main + 10 * np.sqrt(self.cov_main) + self.mean_main / 10,
                                    self.mean_main / 10),
                     color='red',
                     normed=True)
            self.logger.debug(self.model_main_comp[1:])

        plt.hist((self.df[LqCoverage.N_MBASE_COLUMN] / self.df[LqCoverage.QLENGTH_COLUMN]),
                 alpha=0.2,
                 bins=np.arange(0,
                                self.mean_main + 10 * np.sqrt(self.cov_main) + self.mean_main / 10,
                                self.mean_main / 10),
                 color='green',
                 normed=True)
        plt.xlabel('Per read coverage')
        plt.ylabel('Probability density')

        ymin, ymax = plt.gca().get_ylim()
        xmin, xmax = plt.gca().get_xlim()

        plt.legend(bbox_to_anchor=(1,1), loc='upper right', borderaxespad=1)
        if fp:
            plt.savefig(fp, bbox_inches="tight")
        else:
            plt.show()
        plt.close()

    def plot_unmapped_frac_terminal(self, fp=None, adp5_pos=None, adp3_pos=None, x_max=145):
        plt.figure(figsize=(12,5))
        plt.subplot(1,2,1)
        t5l, t3l, il = self.__region_analysis(3, 1)
        self.logger.info("Coordinates of coverage analysis were parsed.")

        plt.hist(t5l, alpha=0.2, bins=np.arange(0, x_max, 5), color='green')
        plt.xlim(0,x_max)
        plt.xlabel('Distance from 5\' terminal')
        plt.ylabel('Frequency')
        ymin5, ymax5 = plt.gca().get_ylim()

        plt.subplot(1,2,2)
        plt.hist(t3l, alpha=0.2, bins=np.arange(0, x_max, 5), color='orange')
        plt.xlim(x_max, 0)
        plt.xlabel('Distance from 3\' terminal')
        plt.ylabel('Frequency')
        ymin3, ymax3 = plt.gca().get_ylim()

        if ymax5 > ymax3:
            ax = plt.subplot(122)
            ax.set_ylim(0, ymax5)
            ymax = ymax5
        else:
            ax = plt.subplot(121)
            ax.set_ylim(0, ymax3)
            ymax = ymax3
        if adp5_pos:
            #adp5_pos -> 45 for pb, 61 for 1d, and 120 for 1d2.
            plt.axvline(x=adp5_pos, linestyle='dashed', linewidth=2, color='red', alpha=0.8) # pacbio
            plt.text(adp5_pos, ymax*0.85, r'Length of the adapter')
        if adp3_pos:
            #adp3_pos -> 45 for pb, 22 for 1d, and 86 for 1d2.
            plt.axvline(x=adp3_pos, linestyle='dashed', linewidth=2, color='red', alpha=0.8) # pacbio
            plt.text(adp3_pos, ymax*0.85, r'Length of the adapter')

        if fp:
            plt.savefig(fp, bbox_inches="tight")
        else:
            plt.show()
        plt.close()

    def plot_qscore_dist(self, platform='ont', fp=None):
        if platform == 'ont':
            mid_threshold = 7 # ont
        else:
            mid_threshold = 8 # pb
        plt.grid(True)
        plt.boxplot([self.df[LqCoverage.QV_COLUMN].values[np.where(self.df[LqCoverage.COVERAGE_COLUMN] == 0.0)],
                     self.df[LqCoverage.QV_COLUMN].values[np.where(self.df[LqCoverage.COVERAGE_COLUMN] != 0.0)]])
        plt.xticks([1,2], ["Non-OLP", "OLP"])
        ymin, ymax = plt.gca().get_ylim()
        plt.axhspan(0,  5, facecolor='red', alpha=0.1)
        plt.axhspan(5,  mid_threshold, facecolor='yellow', alpha=0.1)
        plt.axhspan(mid_threshold, ymax, facecolor='green', alpha=0.1)
        #plt.boxplot(df[5].values[np.where(df[4] == 0.0)])
        plt.ylim(0, ymax)
        plt.ylabel('Averaged QV')
        if fp:
            plt.savefig(fp, bbox_inches="tight")
        else:
            plt.show()
        plt.close()

    def plot_length_vs_coverage(self, fp=None, interval=3000.0):
        ### read score after size binning.
        subplot  = self.__gen_boxplot_length_vs_coverage(interval)
        bin_size = self.df.groupby('Interval').size()
        boundary_reliable_bin = np.where(bin_size >= 50)[0].max()
        xmin, xmax = plt.gca().get_xlim()
        plt.axvspan(boundary_reliable_bin+1.5, xmax+1, facecolor='gray', alpha=0.1)
        #plt.axhline(y=self.mean_main, linestyle='dashed', linewidth=2, color='red', alpha=0.2) # a bit misleading in case skewed dist
        plt.title("Read coverage over different length reads")
        plt.xticks(np.arange(xmax+1), [int(i) for i in np.arange(xmax+1)*interval])
        plt.ylim(0, self.mean_main + 20*np.sqrt(self.cov_main))
        plt.suptitle("")
        #plt.savefig('Box_plot_quality.png')
        if fp:
            plt.savefig(fp, bbox_inches="tight")
        else:
            plt.show()
        plt.close()

        #interval = est_confidence_interval_params(df, 2, nbs=1000)
        #print(interval)
        #print(st.norm.ppf(0.2, loc=interval[0][1], scale=interval[1][0]))

    def __gen_boxplot_length_vs_coverage(self, interval):
        self.df['Interval'] = np.floor(self.df[LqCoverage.QLENGTH_COLUMN].values/interval)
        return self.df.boxplot(column=LqCoverage.COVERAGE_COLUMN, by='Interval', sym='+', rot=90)

    # est_coverage_dist sometimes return weird mu, maybe due to local maxima, so
    #  below method is an ugly working solution.
    # nbs      : the number of bootstrap iteration
    # k_i      : initial value for the parameter for the number of component in GMM.
    # k_max    : maximum value for the parameter for the number of component in GMM.
    def __est_confidence_interval_params(self, nbs=100, k_i=2, k_max=3):
        _mus    = np.zeros(nbs)
        _n      = np.zeros(nbs)
        _sigmas = np.zeros(nbs)

        i = 0
        while(i < nbs):
            (n, mu, sigma_sq) = __est_coverage_dist(k_i=k_i, k_max=k_max)
            _n[i]   = n
            _mus[i] = mu
            _sigmas[i] = np.sqrt(sigma_sq)
            i += 1

        return (st.t.interval(0.99, len(_mus)-1, loc=np.mean(_mus), scale=st.sem(_mus)),
                st.t.interval(0.99, len(_sigmas)-1, loc=np.mean(_sigmas), scale=st.sem(_sigmas)))

    # k_i      : initial value for the parameter for the number of component in GMM.
    # k_max    : maximum value for the parameter for the number of component in GMM.
    def __est_coverage_dist_gmm(self, k_i=2, k_max=10):

        for k in range(k_i, k_max+1):
            _c_i = -1 # index for coverage component. Assume there is only one such a component.
            _b_i = [] # index for background component 

            m_f   = mixture.GaussianMixture(n_components=k).fit(self.df[LqCoverage.COVERAGE_COLUMN].values[np.nonzero(self.df[LqCoverage.COVERAGE_COLUMN])].reshape(-1,1),1)
            order = m_f.weights_/[e[0] for inner in m_f.covariances_ for e in inner] # w/\sigma
            ratio = np.array([e for i in m_f.means_ for e in i])/np.array([e[0] for i in m_f.covariances_ for e in i])
            weird_components = [i for i, v in enumerate(ratio) if v > 1] # sqrt(mu) > sigma

            self.logger.info("The order of componens %s " % " ".join([str(v) for v in order]) )
            self.logger.info("Means of components: %s k=%d" % (" ".join([str(e) for i in m_f.means_ for e in i]), k))
            self.logger.info("Covariances of components: %s k=%d" % (" ".join([str(e[0]) for inner in m_f.covariances_ for e in inner]), k))

            _max = -1*np.inf
            for i,v in enumerate(order):
                if i in weird_components:
                    continue

                if _max < v:
                    _max = v
                    _c_i = i

            for i in range(0, len(order)):
                if i in weird_components or i == _c_i:
                    continue
                else:
                    _b_i.append(i)
        
            if _c_i == -1 or not _b_i:
                continue
            else:
                return (m_f, m_f.means_[_c_i][0], m_f.covariances_[_c_i][0][0])
                #return (k, m_f, m_f.means_[_c_i][0], m_f.covariances_[_c_i][0][0])

    def __region_analysis(self, column_i_coords, column_i_ql, threshold=50):

        coi = column_i_coords
        qli = column_i_ql

        trim_5 = []
        trim_3 = []
        intrnl = []

        for i, str in enumerate(self.df[coi]):
            if str == '0':
                continue
            ql = self.df[qli][i]
            regs = [ (int(reg.split('-')[0]), int(reg.split('-')[1])) for reg in str.split(',')]
            s = e = None
            if(len(regs) > 1):
                sr = sorted(regs, key=itemgetter(0,1))
                s = sr[0][0]
                e = sr[-1][1]
                for k in range(0, len(regs)-1):
                    intrnl.append(regs[k+1][0] - regs[k][1])
            elif(len(regs) == 1):
                (s, e) = regs[0]
            else:
                self.logger.warning("The number of region is weird." )

            if s != None and e != None:
                t3 = int(ql) - int(e)
                trim_5.append(s)
                trim_3.append(t3)

        return(trim_5, trim_3, intrnl)


# test
if __name__ == "__main__":

    lc = LqCoverage("/home/fukasay/temp/ont_ecoli_k12w5m40_oh2000_r40_min1000_debug.txt")
    lc.plot_coverage_dist()
    lc.plot_unmapped_frac_terminal(adp5_pos=61, adp3_pos=None)
    lc.plot_qscore_dist()
    lc.plot_length_vs_coverage()

    """
    # internal break. experimental.
    plt.hist(il, alpha=0.2, color='blue')
    print(len(il))
    plt.plot(il, ([-0.5]*len(il)) + 0.1*np.random.random(len(il)), '+k')
    plt.show()
    """

