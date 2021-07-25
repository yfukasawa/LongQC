import math, logging, sys
import numpy             as np
import pandas            as pd
import scipy.stats       as st
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from sklearn  import mixture
from operator import itemgetter
from scipy.signal import argrelmax
from mixEM import mixem

from logging import getLogger
logger = getLogger(__name__)

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
    COV_CORRECTION              = 0.9 # experimental term.
    DIV_SCORE_THRESHOLD         = 0.25 # experimental.
    COV_THRESHOLD_FOR_DIV_SC    = 25  # experimental.
    LENGTH_BIN_THRESHOLD        = 100
    # table index
    READ_NAME_COLUMN   = 0
    QLENGTH_COLUMN     = 1
    N_MBASE_COLUMN     = 2
    MED_READ_COV_CORS  = 4
    #GOOD_READ_COV_CORS = 5
    T1_COVERAGE_COLUMN = 5
    QV_COLUMN          = 6
    DIV_COLUMN         = 7
    COVERAGE_COLUMN    = 8

    def __init__(self, table_path, isTranscript=False, control_filtering=None, engine='python'):
        self.df = pd.read_table(table_path, sep='\t', header=None, dtype={3: str, 4: str})
        self.min_lambda = None
        self.max_lambda = None
        self.unmapped_frac_trimmed   = -1.0
        self.unmapped_frac_untrimmed = -1.0
        self.unmapped_frac_med       = -1.0
        #self.unmapped_bad_frac       = -1.0
        self.high_div_frac           = -1.0
        self.model = None
        self.mean_main = None
        self.cov_main  = None
        self.main_comp_index = None
        self.control_reads = None
        self.low_coverage = None
        self.no_coverage = None

        if control_filtering is not None:
            self.df_control = pd.read_table(control_filtering, sep='\t', header=None)
            self.control_reads = self.df_control[self.df_control[LqCoverage.T1_COVERAGE_COLUMN] >= 0.5][0].tolist()
            self.df = self.df[~self.df[LqCoverage.READ_NAME_COLUMN].isin(self.control_reads)] # remove control from the list

        # for final conclusion
        self.warnings = []
        self.errors = []

        # for transcript type data. e.g. direct RNA, cDNA, Iso-seq, etc
        self.mix_model = None
        self.mode_logn_main  = None
        self.mu_logn_main    = None
        self.sigma_logn_main = None
        
        self.isTranscript = isTranscript
        
        logger.info("Estimating coverage distribution..")
        self.__est_coverage()
        logger.info("Estimation of coverage distribution finished.")

    def get_mean(self):
        if not self.mean_main:
            logger.warning("Mean has no value. Do estimation first.")
        return self.mean_main

    def get_sd(self):
        if not self.cov_main:
            logger.warning("SD has no value. Do estimation first.")
        return np.sqrt(self.cov_main)

    def get_logn_mode(self):
        if not self.mode_logn_main:
            logger.warning("Mode of lognormal has no value. Do estimation first.")
            return None
        return self.mode_logn_main

    def get_logn_mu(self):
        if not self.mu_logn_main:
            logger.warning("Mu of lognormal has no value. Do estimation first.")
        return self.mu_logn_main

    def get_logn_sigma(self):
        if not self.sigma_logn_main:
            logger.warning("sigma of lognormal has no value. Do estimation first.")
        return self.sigma_logn_main

    def get_expected_zero_rate(self):
        if not self.mode_logn_main and not self.mean_main:
            logger.warning("Mode of lognormal and mean of GMM has no value. Do estimation first.")
            return None
        elif not self.mode_logn_main:
            return (self.mean_main, 1.3865*0.64086**self.mean_main)
        else:
            return (self.mode_logn_main, 1.3865*0.64086**self.mode_logn_main) # experimental

    def get_unmapped_frac(self):
        if self.unmapped_frac_trimmed == -1.0:
            logger.warning("Unmapped fraction has no value. Do estimation first.")
        return self.unmapped_frac_trimmed

    def get_unmapped_med_frac(self):
        if self.unmapped_frac_med == -1.0:
            logger.warning("Unmapped medial fraction has no value. Do estimation first.")
        return self.unmapped_frac_med

    #def get_unmapped_bad_frac(self):
    #    if self.unmapped_bad_frac == -1.0:
    #        logger.warning("Unmapped bad read fraction has no value. Do estimation first.")
    #    return self.unmapped_bad_frac

    def is_no_coverage(self):
        if self.no_coverage != None:
            return self.no_coverage
        else:
            return None

    def is_low_coverage(self):
        if self.low_coverage != None:
            return self.low_coverage
        else:
            return None

    # experimental
    def get_high_div_frac(self):
        if self.high_div_frac == -1.0:
            logger.warning("Highly divergent read fraction has no value. Do estimation first.")
        return self.high_div_frac

    def get_control_num(self):
        if self.control_reads:
            return len(self.control_reads)
        else:
            return 0.0

    def get_control_frac(self):
        if self.control_reads:
            return len(self.control_reads)/(len(self.control_reads) + len(self.df))
        else:
            return 0.0

    def get_errors(self):
        return self.errors

    def get_warnings(self):
        return self.warnings

    def __est_coverage(self):
        self.unmapped_frac_trimmed   = self.df[LqCoverage.T1_COVERAGE_COLUMN].values[np.where(self.df[LqCoverage.T1_COVERAGE_COLUMN] == 0.0)].shape[0] \
                                       / self.df.shape[0]
        self.unmapped_frac_untrimmed = self.df[LqCoverage.N_MBASE_COLUMN].values[np.where(self.df[LqCoverage.N_MBASE_COLUMN] == 0.0)].shape[0] \
                                       / self.df.shape[0]
        self.unmapped_frac_med       = self.df[LqCoverage.MED_READ_COV_CORS].values[np.where(self.df[LqCoverage.MED_READ_COV_CORS] == '0')].shape[0] \
                                       / self.df.shape[0]
        #self.unmapped_bad_frac       = self.df[LqCoverage.GOOD_READ_COV_CORS].values[np.where(self.df[LqCoverage.GOOD_READ_COV_CORS] == '0')].shape[0] \
        #                               / self.df.shape[0]
        self.high_div_frac           = self.df[LqCoverage.DIV_COLUMN].values[\
                                       np.where((self.df[LqCoverage.DIV_COLUMN] >= LqCoverage.DIV_SCORE_THRESHOLD) & \
                                                (self.df[LqCoverage.T1_COVERAGE_COLUMN] >= LqCoverage.COV_THRESHOLD_FOR_DIV_SC) & \
                                                (self.df[LqCoverage.MED_READ_COV_CORS] != '0') )].shape[0] \
                                       / self.df.shape[0]
        #logger.info("Unmapped fraction: %.3f (naive), %.3f (coverage considered)" % (self.unmapped_frac_untrimmed, self.unmapped_frac_trimmed))

        model_main_comp = self.__est_coverage_dist_gmm(k_i=2)
            
        self.model = model_main_comp[0]
        self.mean_main = model_main_comp[1]
        self.cov_main  = model_main_comp[2]
        self.main_comp_index = model_main_comp[3]

        raw_hist = plt.hist((self.df[LqCoverage.N_MBASE_COLUMN] / self.df[LqCoverage.QLENGTH_COLUMN]),
                            alpha=0.2,
                            bins=np.arange(0,
                                           self.mean_main + 10 * np.sqrt(self.cov_main) + self.mean_main / 10,
                                           self.mean_main / 10),
                            color='green',
                            density=True)
        plt.close()

        self.low_coverage = self.__looks_lowcoverage(raw_hist)

        if self.unmapped_frac_med >= LqCoverage.UNMAPPED_FRACTION_THRESHOLD:
            logger.warning("The fraction of zero coverage read is high %.3f" % self.unmapped_frac_med)
            self.min_lambda = -1 * math.log(self.unmapped_frac_med - LqCoverage.UNMAPPED_FRACTION_PARAM_MIN) 
            self.max_lambda = -1 * math.log(self.unmapped_frac_med - LqCoverage.UNMAPPED_FRACTION_PARAM_MAX)
            range_str = str(self.min_lambda) + "-" + str(self.max_lambda)
            logger.warning("If and only if the data is healthy, very rough estimated coverage range is %s." % range_str)

        if self.model is None:
            self.low_coverage = None # no_coverage is a stronger condition.
            self.no_coverage = True
            logger.warning("No coverage data is available.")
            return

        # skewed dist. sometimes occur in low coverage data and/or long tail
        if self.low_coverage and self.unmapped_frac_med < LqCoverage.UNMAPPED_FRACTION_THRESHOLD and not self.isTranscript:
            logger.warning("The fraction of zero coverage read is not too high %.3f, but still looks low coverage" % self.unmapped_frac_med)
            logger.warning(self.unmapped_frac_med)

            self.__est_coverage_dist_lognorm_norm()
            self.mode_logn_main = np.exp(self.mix_model[1][1] - self.mix_model[2][1]**2)
            self.mu_logn_main = self.mix_model[1][1]
            self.sigma_logn_main = self.mix_model[2][1]

	# edgy case fix.
        if self.low_coverage and self.unmapped_frac_med >= LqCoverage.UNMAPPED_FRACTION_THRESHOLD and not self.isTranscript:
            logger.warning("The fraction of zero coverage read is too high %.3f." % self.unmapped_frac_med)
            logger.warning(self.unmapped_frac_med)

            self.__est_coverage_dist_lognorm_norm()
            self.mode_logn_main = np.exp(self.mix_model[1][1] - self.mix_model[2][1]**2)
            self.mu_logn_main = self.mix_model[1][1]
            self.sigma_logn_main = self.mix_model[2][1]

        # RNA and long tail case
        if self.isTranscript:
            logger.warning("Distribution looks long tail RNA type coverage. Apply lognorm.")
            self.__est_coverage_dist_lognorm_norm()
            self.mode_logn_main = np.exp(self.mix_model[1][1] - self.mix_model[2][1]**2*0.5)
            self.mu_logn_main = self.mix_model[1][1]
            self.sigma_logn_main = self.mix_model[2][1]


    def __looks_lowcoverage(self, hist):
        relmaxs = argrelmax(hist[0])
        if hist[0][0] / np.sum(hist[0]) <0.01:
            return False
        for mx_i in relmaxs[0]:
            if hist[0][mx_i] > (hist[0][0] / 5):
                return False
                break
        return True

    def plot_coverage_dist(self, fp=None):
        if self.min_lambda and self.max_lambda:
            plt.figure(figsize=(12,5))
            ax1 = plt.subplot(1,2,1)
            plt.grid(True)
            pois_min = st.poisson(self.min_lambda)
            pois_max = st.poisson(self.max_lambda)
            pois_x = np.arange(int(self.mean_main+4 * np.sqrt(self.cov_main)) + 1)
            pois_y_min = pois_min.pmf(pois_x)
            pois_y_max = pois_max.pmf(pois_x)
            plt.xlim(0, 50)
            plt.plot(pois_x, pois_y_min, label="Fitted Model by Poisson model (" + "%.3f" % self.min_lambda + ")")
            plt.plot(pois_x, pois_y_max, label="Fitted Model by Poisson model (" + "%.3f" % self.max_lambda + ")")
            plt.xlabel('Per read coverage')
            plt.ylabel('Probability density')
            plt.hist(self.df[LqCoverage.COVERAGE_COLUMN],
                 alpha=0.2,
                 bins=np.arange(0,
                                self.mean_main + 10 * np.sqrt(self.cov_main) + self.mean_main / 10,
                                self.mean_main / 10),
                 color='green',
                 density=True)
            plt.legend(bbox_to_anchor=(1,1), loc='upper right', borderaxespad=1)
            plt.subplot(1,2,2)
            plt.grid(True)
        else:
            plt.grid(True)

        gmm_x = np.linspace(0, self.mean_main+10*np.sqrt(self.cov_main), 5000)
        if self.mix_model is not None:
            # RNA case. 0 -> bg and 1 -> target
            mix_y = self.mix_model[0][0] * st.norm(self.mix_model[1][0], self.mix_model[2][0]).pdf(gmm_x) + \
                    self.mix_model[0][1] * st.lognorm.pdf(gmm_x, self.mix_model[2][1], loc=0, scale=np.exp(self.mix_model[1][1]))
            plt.plot(gmm_x, mix_y, label="Fitted by Lognormal and gaussian mixture model")
            plt.xlim(0, self.mean_main+10*np.sqrt(self.cov_main))
            plt.legend(bbox_to_anchor=(1,1), loc='upper right', borderaxespad=1)
        elif self.model is None:
            plt.gcf().text(0.30,0.5,"Caution: coverage estimation was skipped due to insufficient amount of data.", backgroundcolor='yellow')            
        else:
            gmm_y = np.exp(self.model.score_samples(gmm_x.reshape(-1,1)))
            plt.plot(gmm_x, gmm_y, label="Fitted by Gaussian mixture model")
            plt.xlim(0, self.mean_main+10*np.sqrt(self.cov_main))
            plt.legend(bbox_to_anchor=(1,1), loc='upper right', borderaxespad=1)

        #plt.hist((self.df[LqCoverage.N_MBASE_COLUMN] / self.df[LqCoverage.QLENGTH_COLUMN]),
        #         alpha=0.2,
        #         bins=np.arange(0,
        #                        self.mean_main + 10 * np.sqrt(self.cov_main) + self.mean_main / 10,
        #                        self.mean_main / 10),
        #         color='green',
        #         normed=True)
        plt.hist(self.df[LqCoverage.COVERAGE_COLUMN],
                 alpha=0.2,
                 bins=np.arange(0,
                                self.mean_main + 10 * np.sqrt(self.cov_main) + self.mean_main / 10,
                                self.mean_main / 10),
                 color='green',
                 density=True)

        plt.xlabel('Per read coverage')
        plt.ylabel('Probability density')

        ymin, ymax = plt.gca().get_ylim()
        xmin, xmax = plt.gca().get_xlim()

        if fp:
            plt.savefig(fp, bbox_inches="tight")
        else:
            plt.show()
        plt.close()

    def calc_xome_size(self, throughput):
        m_size = -1
        if self.no_coverage:
            return "N/A"
        elif self.isTranscript:
            m_size = int((throughput * (1.0 - self.unmapped_frac_med)) / self.mode_logn_main)
        elif self.low_coverage:
            m_size = int((throughput * (1.0 - self.unmapped_frac_med)) / self.mode_logn_main)
        else:
            # gmm
            m_size = int((throughput * (1.0 - self.unmapped_frac_med)) / self.mean_main)

        if self.unmapped_frac_med >= LqCoverage.UNMAPPED_FRACTION_THRESHOLD: #or (self.low_coverage and not self.isTranscript):
            # poission est + empirical correction
            _s1 = throughput * LqCoverage.COV_CORRECTION * (1 - LqCoverage.UNMAPPED_FRACTION_PARAM_MIN) / self.min_lambda
            _s2 = throughput * LqCoverage.COV_CORRECTION * (1 - LqCoverage.UNMAPPED_FRACTION_PARAM_MAX) / self.max_lambda
            return "%d (e = %.1f%%), %d (e = 20%%), %d (e = 5%%)" % (m_size, self.unmapped_frac_med*100, _s2, _s1)
        else:
            return "%d (e = %.1f%%)" % (m_size, self.unmapped_frac_med*100)

    def plot_unmapped_frac_terminal(self, fp=None, *, adp5_pos=None, adp3_pos=None, x_max=145):
        plt.figure(figsize=(12,5))
        ax1 = plt.subplot(1,2,1)
        t5l, t3l, il = self.__region_analysis(3, 1)
        logger.info("Coordinates of coverage analysis were parsed.")

        plt.axes(ax1)
        plt.hist(t5l, alpha=0.2, bins=np.arange(0, x_max, 5), color='green')
        plt.xlim(0,x_max)
        plt.xlabel('Distance from 5\' terminal')
        plt.ylabel('Frequency')
        ymin5, ymax5 = plt.gca().get_ylim()

        ax2 = plt.subplot(1,2,2)
        plt.axes(ax2)
        plt.hist(t3l, alpha=0.2, bins=np.arange(0, x_max, 5), color='orange')
        plt.xlim(x_max, 0)
        plt.xlabel('Distance from 3\' terminal')
        plt.ylabel('Frequency')
        ymin3, ymax3 = plt.gca().get_ylim()

        if ymax5 > ymax3:
            ax2.set_ylim(0, ymax5)
            ymax = ymax5
        else:
            ax1.set_ylim(0, ymax3)
            ymax = ymax3

        if adp5_pos:
            #adp5_pos -> 45 for pb, 61 for 1d, and 120 for 1d2.
            ax1.axvline(x=adp5_pos, linestyle='dashed', linewidth=2, color='red', alpha=0.8) # pacbio
            if adp5_pos > 90:
                ax1.text(adp5_pos, ymax*0.85, r'Length of the adapter', horizontalalignment='right')
            else:
                ax1.text(adp5_pos, ymax*0.85, r'Length of the adapter', horizontalalignment='left')

        if adp3_pos:
            #adp3_pos -> 45 for pb, 22 for 1d, and 86 for 1d2.
            ax2.axvline(x=adp3_pos, linestyle='dashed', linewidth=2, color='red', alpha=0.8) # pacbio
            if adp3_pos > 90:
                ax2.text(adp3_pos, ymax*0.85, r'Length of the adapter', horizontalalignment='left')
            else:
                ax2.text(adp3_pos, ymax*0.85, r'Length of the adapter', horizontalalignment='right')

        if fp:
            plt.savefig(fp, bbox_inches="tight")
        else:
            plt.show()
        plt.close()

    def plot_qscore_dist(self, fp=None, *, platform='ont'):
        if platform == 'ont':
            mid_threshold = 7 # ont
        else:
            #id_threshold = 8 # obsolete
            mid_threshold = 7 # we can use the same criterion for pb too.
        plt.grid(True)
        plt.boxplot([self.df[LqCoverage.QV_COLUMN].values[np.where(self.df[LqCoverage.COVERAGE_COLUMN] == 0.0)],
                     self.df[LqCoverage.QV_COLUMN].values[np.where(self.df[LqCoverage.COVERAGE_COLUMN] != 0.0)]])
        plt.xticks([1,2], ["Non-sense reads", "Normal reads"])
        ymin, ymax = plt.gca().get_ylim()
        #plt.axhspan(0,  5, facecolor='red', alpha=0.1)
        #plt.axhspan(5,  mid_threshold, facecolor='yellow', alpha=0.1)
        plt.axhspan(0, mid_threshold, facecolor='red', alpha=0.1)
        plt.axhspan(mid_threshold, ymax, facecolor='green', alpha=0.1)
        #plt.boxplot(df[5].values[np.where(df[4] == 0.0)])
        plt.ylim(0, ymax)
        plt.ylabel('Averaged QV')
        if fp:
            plt.savefig(fp, bbox_inches="tight")
        else:
            plt.show()
        plt.close()

    def plot_length_vs_coverage(self, fp=None, *, interval=3000.0):
        ### read score after size binning.
        subplot  = self.__gen_boxplot_length_vs_coverage(interval)
        bin_size = self.df.groupby('Binned read length').size()
        boundary_reliable_bins = np.where(bin_size >=  LqCoverage.LENGTH_BIN_THRESHOLD)[0]
        xmin, xmax = plt.gca().get_xlim()
        if boundary_reliable_bins.size > 0:
            dmin = boundary_reliable_bins.min()
            dmax = boundary_reliable_bins.max()
            if dmax < xmax:
                plt.axvspan(boundary_reliable_bins.max()+1.5, xmax+1, facecolor='gray', alpha=0.1)
            if dmin > xmin:
                plt.axvspan(xmin-1, dmin+1.5, facecolor='gray', alpha=0.1)
        else:
            plt.axvspan(xmin-1,xmax+1,facecolor='gray', alpha=0.1)
        plt.xlim(xmin, xmax)
        #plt.axhline(y=self.mean_main, linestyle='dashed', linewidth=2, color='red', alpha=0.2) # a bit misleading in case skewed dist
        plt.title("Read coverage over different length reads")
        plt.xticks(np.arange(xmax+1), [int(i) for i in np.arange(xmax+1)*interval])
        plt.ylim(0, self.mean_main + 20*np.sqrt(self.cov_main))
        plt.ylabel("per-read coverage")
        plt.suptitle("")

        if not self.min_lambda and not self.max_lambda and self.mix_model is None:
            # 3 sigma
            yc = self.get_mean() - 3*self.get_sd()
            plt.axhline(y=yc, color="royalblue", alpha=0.4, lw=1)
            plt.text(0, yc, r'3$\sigma$', color="royalblue")
            yc = self.get_mean() + 3*self.get_sd()
            plt.axhline(y=yc, color="royalblue", alpha=0.4, lw=1)
            plt.text(0, yc, r'3$\sigma$', color="royalblue")

            self.__check_outlier_coverage(interval)

        #plt.savefig('Box_plot_quality.png')
        if fp:
            plt.savefig(fp, bbox_inches="tight", transparent=True)
        else:
            plt.show()
        plt.close()
        #interval = est_confidence_interval_params(df, 2, nbs=1000)
        #print(interval)
        #print(st.norm.ppf(0.2, loc=interval[0][1], scale=interval[1][0]))

    def __gen_boxplot_length_vs_coverage(self, interval):
        self.df.loc[self.df[LqCoverage.QLENGTH_COLUMN] >= 3000, 'MERGED_COVERAGE'] = self.df[LqCoverage.COVERAGE_COLUMN]
        self.df.loc[self.df[LqCoverage.QLENGTH_COLUMN] <  3000, 'MERGED_COVERAGE'] = self.df[LqCoverage.T1_COVERAGE_COLUMN]
        self.df['Binned read length'] = np.floor(self.df[LqCoverage.QLENGTH_COLUMN].values/interval)
        #return self.df.boxplot(column=LqCoverage.COVERAGE_COLUMN, by='Binned read length', sym='+', rot=90, figsize=(2*int(max(self.df['Binned read length'])/5+0.5), 4.8))
        #return self.df.boxplot(column='MERGED_COVERAGE', by='Binned read length', sym='+', rot=90, figsize=(2*int(max(self.df['Binned read length'])/5+0.5), 4.8))
        if max(self.df['Binned read length']) < 5:
            return self.df.boxplot(column='MERGED_COVERAGE', by='Binned read length', sym='+', rot=90)
        else:
            return self.df.boxplot(column='MERGED_COVERAGE', by='Binned read length', sym='+', rot=90, figsize=(2*int(max(self.df['Binned read length'])/5+0.5), 4.8))

    def __check_outlier_coverage(self, interval):
        stats = self.df.groupby('Binned read length')[LqCoverage.COVERAGE_COLUMN].agg([np.median, np.size])
        meds = stats['median'].iloc[np.where(stats['size']>=LqCoverage.LENGTH_BIN_THRESHOLD)[0]]
        three_sigma = np.where((meds > self.get_mean() + 3*self.get_sd()) | (meds <= self.get_mean() - 3*self.get_sd()))
        if len(three_sigma[0]) > 0:
            #error case
            #self.errors.append(('Coverage error', 'Coverage is not homogenous over the read length.'))
            # change this to warnings
            self.warnings.append(('Coverage warning', 'Coverage might not be homogenous over the read length.'))
        #else:
        #    two_sigma = np.where((meds > self.get_mean() + 2*self.get_sd()) | (meds <= self.get_mean() - 2*self.get_sd()))
        #    if len(two_sigma[0]):
        #        self.warnings.append(('Coverage warning', 'Coverage might not be homogenous over the read length.'))

    # est_coverage_dist sometimes return weird mu, maybe due to local maxima, so
    #  below method is an ugly working solution.
    # nbs      : the number of bootstrap iteration
    # k_i      : initial value for the parameter for the number of component in GMM.
    # k_max    : maximum value for the parameter for the number of component in GMM.
    def __est_confidence_interval_params(self, nbs=100, k_i=2, k_max=2):
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

    def __est_coverage_dist_lognorm_norm(self):
        #a little bit buggy
        th_per    = self.df[LqCoverage.COVERAGE_COLUMN].quantile(0.85)
        if th_per == 0.0:
            th_per    = self.df[LqCoverage.COVERAGE_COLUMN].quantile(1.0)
        nonzeros  = self.df[LqCoverage.COVERAGE_COLUMN].values[self.df[LqCoverage.COVERAGE_COLUMN].to_numpy().nonzero()]
        # we should not assume k=2 model, but tentatively work with this.
        i_bg = 0 if self.main_comp_index == 1 else 1
        i_m  = 1 if self.main_comp_index == 1 else 0
        weights, distributions, ll = mixem.em(nonzeros[nonzeros < th_per], [
            mixem.distribution.NormalDistribution(self.model.means_[i_bg][0], np.sqrt(self.model.covariances_[i_bg][0][0])), # for noise
            mixem.distribution.LogNormalDistribution(np.log(self.model.means_[i_m][0]), 1) # for true
        ], max_iterations=500)

        self.mix_model = (weights, [d.get_mu() for d in distributions], [d.get_sigma() for d in distributions])

    # k_i      : initial value for the parameter for the number of component in GMM.
    # k_max    : maximum value for the parameter for the number of component in GMM.
    def __est_coverage_dist_gmm(self, k_i=2, k_max=2):

        _c_i = -1 # index for coverage component. Assume there is only one such a component.
        _b_i = [] # index for background component 
        for k in range(k_i, k_max+1):
            _c_i = -1 # init
            _b_i = [] # init

            #a little bit buggy
            th_per = self.df[LqCoverage.COVERAGE_COLUMN].quantile(0.85)
            if th_per == 0.0:
                th_per    = self.df[LqCoverage.COVERAGE_COLUMN].quantile(1.0)
            nonzeros  = self.df[LqCoverage.COVERAGE_COLUMN].values[self.df[LqCoverage.COVERAGE_COLUMN].to_numpy().nonzero()]
            if nonzeros[nonzeros < th_per].size == 0:
                # literally no data case!
                # give some dummy
                return (None, 1, 10, 0)
            else:
                m_f   = mixture.GaussianMixture(n_components=k).fit(nonzeros[nonzeros < th_per].reshape(-1,1),1)

            logger.debug(m_f)
            order = m_f.weights_/[e[0] for inner in m_f.covariances_ for e in inner] # w/\sigma
            ratio = np.array([e for i in m_f.means_ for e in i])/np.array([e[0] for i in m_f.covariances_ for e in i])
            #weird_components = [i for i, v in enumerate(ratio) if v > 1] # sqrt(mu) > sigma

            logger.info("The order of componens %s " % " ".join([str(v) for v in order]) )
            logger.info("Means of components: %s k=%d" % (" ".join([str(e) for i in m_f.means_ for e in i]), k))
            logger.info("Covariances of components: %s k=%d" % (" ".join([str(e[0]) for inner in m_f.covariances_ for e in inner]), k))

            _max = -1*np.inf
            for i,v in enumerate(order):
                #if i in weird_components:
                #    continue

                if _max < v:
                    _max = v
                    _c_i = i

            for i in range(0, len(order)):
                #if i in weird_components or i == _c_i:
                if i == _c_i:
                    continue
                else:
                    _b_i.append(i)
        
            if _c_i == -1 or not _b_i:
                continue
            else:
                return (m_f, m_f.means_[_c_i][0], m_f.covariances_[_c_i][0][0], _c_i)
                #return (k, m_f, m_f.means_[_c_i][0], m_f.covariances_[_c_i][0][0])

        return (m_f, m_f.means_[_c_i][0], m_f.covariances_[_c_i][0][0], _c_i)

    def __region_analysis(self, column_i_coords, column_i_ql, threshold=50):

        coi = column_i_coords
        qli = column_i_ql

        trim_5 = []
        trim_3 = []
        intrnl = []

        for i in self.df.index.tolist():
            str = self.df[coi][i]
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
                logger.warning("The number of region is weird." )

            if s != None and e != None:
                t3 = int(ql) - int(e)
                trim_5.append(s)
                trim_3.append(t3)

        return(trim_5, trim_3, intrnl)

# test
if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(
        prog='lq_coverage.py',
        description='a module for LongQC. gestimates coverage and generates plots and some stats from tables.',
        add_help=True,
    )
    parser.add_argument('input', \
                        help='path for coverage_out.txt made by minimap2-coverage.', type=str)
    parser.add_argument('--output', dest='outf', \
                        help='output path where plots are saved.', type=str)
    parser.add_argument('--control', dest='conf',\
                        help='(optional) path for spikein_out.txt.', type=str)
    args = parser.parse_args()

    if args.conf:
        conf = args.conf
    else:
        conf = None

    inf  = args.input
    outf = args.outf

    lc = LqCoverage(inf, isTranscript=False, control_filtering=conf)
    if lc.mode_logn_main:
        print("Mode of LogN mix: %.3f" % lc.mode_logn_main)
        print("Mean of GMM: %.3f" % lc.mean_main)
    else:
        print("Mean of GMM: %.3f" % lc.mean_main)
    plt.rcParams['figure.figsize'] = (7, 7)
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42

    # comment out a plot if needed. not multiple plots in a single file. 
    #lc.plot_qscore_dist(outf)
    lc.plot_coverage_dist(outf)
    #lc.plot_length_vs_coverage(outf)
    #lc.plot_unmapped_frac_terminal(outf)

    # logs of stats.
    print("%% non-sense reads: %.3f" % lc.get_unmapped_med_frac())
    print("%% control reads: %.3f" % lc.get_control_frac())
    print("Warnings: ", lc.warnings)
    print("Errors: ", lc.errors)

    """
    # internal break. experimental.
    plt.hist(il, alpha=0.2, color='blue')
    print(len(il))
    plt.plot(il, ([-0.5]*len(il)) + 0.1*np.random.random(len(il)), '+k')
    plt.show()
    """

