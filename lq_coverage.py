import math
import numpy             as np
import pandas            as pd
import scipy.stats       as st
import matplotlib.pyplot as plt

#from scipy.stats   import gumbel_r
from sklearn  import mixture
from lq_utils import eprint
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
"""
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


def gen_boxplot_length_vs_coverage(df, interval):
    df['Interval'] = np.floor(df[1].values/interval)
    return df.boxplot(column=4, by='Interval', sym='+', rot=90)


# est_coverage_dist sometimes return weird mu, maybe due to local maxima, so
#  below method is an ugly working solution.

# df       : DataFrame of pandas
# column_i : index for coverage in df
# nbs      : the number of bootstrap iteration
# k_i      : initial value for the parameter for the number of component in GMM.
# k_max    : maximum value for the parameter for the number of component in GMM.
def est_confidence_interval_params(df, column_i, nbs=100, k_i=2, k_max=3):
    _mus    = np.zeros(nbs)
    _n      = np.zeros(nbs)
    _sigmas = np.zeros(nbs)

    i = 0
    while(i < nbs):
        (n, mu, sigma_sq) = est_coverage_dist(df, column_i, k_i=k_i, k_max=k_max)
        _n[i]   = n
        _mus[i] = mu
        _sigmas[i] = np.sqrt(sigma_sq)
        i += 1

    return (st.t.interval(0.99, len(_mus)-1, loc=np.mean(_mus), scale=st.sem(_mus)),
            st.t.interval(0.99, len(_sigmas)-1, loc=np.mean(_sigmas), scale=st.sem(_sigmas)))


# df       : DataFrame of pandas
# column_i : index for coverage in df
# nbs      : the number of bootstrap iteration
# k_i      : initial value for the parameter for the number of component in GMM.
# k_max    : maximum value for the parameter for the number of component in GMM.
def est_coverage_dist(df, column_i, k_i=2, k_max=10):

    for k in range(k_i, k_max+1):
        _c_i = -1 # index for coverage component. Assume there is only one such a component.
        _b_i = [] # index for background component 

        m_f   = mixture.GaussianMixture(n_components=k).fit(df[column_i].values[np.nonzero(df[column_i])].reshape(-1,1),1)
        #m_f   = mixture.GaussianMixture(n_components=k).fit(df[column_i].values.reshape(-1,1),1)
        order = m_f.weights_/[e[0] for inner in m_f.covariances_ for e in inner] # w/\sigma
        ratio = np.array([e for i in m_f.means_ for e in i])/np.array([e[0] for i in m_f.covariances_ for e in i])
        weird_components = [i for i, v in enumerate(ratio) if v > 1] # sqrt(mu) > sigma

        #print("debug", weird_components)
        print("debug", order)
        print("debug", [e for i in m_f.means_ for e in i], "k =", k)
        print("debug", [e[0] for inner in m_f.covariances_ for e in inner], "k =", k)

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


def region_analysis(df, column_i_coords, column_i_ql, threshold=50):

    coi = column_i_coords
    qli = column_i_ql

    trim_5 = []
    trim_3 = []
    intrnl = []

    for i, str in enumerate(df[coi]):
        if str == '0':
            continue
        ql = df[qli][i]
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
            eprint("Error: the number of region is weird." )

        if s != None and e != None:
            t3 = int(ql) - int(e)
            trim_5.append(s)
            trim_3.append(t3)

    return(trim_5, trim_3, intrnl)


# test
if __name__ == "__main__":
    df = pd.read_table('/home/fukasay/temp/sequel_arabi_k12w5m40_ava_sub10k_coverage_minimap2mod_bothfilters_oh2000_ratio40_debug.txt', sep='\t', header=None)

    unmapped_fraction_threshold = 0.4
    unmap_fraction_param_min    = 0.05
    unmap_fraction_param_max    = 0.2

    coverage_column_i = 4

    unmapped_frac_trimmed   = len(df[coverage_column_i].values[np.where(df[coverage_column_i] == 0.0)])/len(df[coverage_column_i])
    unmapped_frac_untrimmed = len(df[2].values[np.where(df[2] == 0.0)])/len(df[1])
    print(unmapped_frac_untrimmed, unmapped_frac_trimmed)
    #print(unmapped_frac_trimmed)

    min_lambda = max_lambda = None
    #mean_main = cov_main = None

    if unmapped_frac_trimmed >= unmapped_fraction_threshold:
        print("Warning: the fraction of zero coverage read is high,", unmapped_frac_trimmed)
        min_lambda = -1 * math.log(unmapped_frac_trimmed-unmap_fraction_param_min)
        max_lambda = -1 * math.log(unmapped_frac_trimmed-unmap_fraction_param_max)
        print("If and only if the data is healthy, very rough estimated coverage range is",str(min_lambda), "-", str(max_lambda)+".")

    model_main_comp = est_coverage_dist(df, coverage_column_i, k_i=2)
    model = model_main_comp[0]
    mean_main = model_main_comp[1]
    cov_main = model_main_comp[2]

    gmm_x = np.linspace(0, mean_main+4*np.sqrt(cov_main), 5000)
    gmm_y = np.exp(model.score_samples(gmm_x.reshape(-1,1)))
    plt.grid(True)

    if min_lambda and max_lambda:
        pois_min = st.poisson(min_lambda)
        pois_max = st.poisson(max_lambda)
        pois_x = np.arange(int(mean_main+4*np.sqrt(cov_main))+1)
        pois_y_min = pois_min.pmf(pois_x)
        pois_y_max = pois_max.pmf(pois_x)
        plt.xlim(0, 50)
        plt.plot(pois_x, pois_y_min, label="Fitted Model by Poisson model (" + "%.3f" % min_lambda + ")")
        plt.plot(pois_x, pois_y_max, label="Fitted Model by Poisson model (" + "%.3f" % max_lambda + ")")
    else:
        plt.plot(gmm_x, gmm_y, label="Fitted by Gaussian mixture model")

    if unmapped_frac_trimmed < unmapped_fraction_threshold:
        plt.hist(df[coverage_column_i], alpha=0.2, bins=np.arange(0, mean_main+10*np.sqrt(cov_main) + mean_main/10, mean_main/10), color='red', normed=True)
        print(model_main_comp[1:])

    plt.hist((df[2]/df[1]), alpha=0.2, bins=np.arange(0, mean_main+10*np.sqrt(cov_main) + mean_main/10, mean_main/10), color='green', normed=True)
    plt.xlabel('Per read coverage')
    plt.ylabel('Probability density')

    ymin, ymax = plt.gca().get_ylim()
    xmin, xmax = plt.gca().get_xlim()

    plt.legend(bbox_to_anchor=(1,1), loc='upper right', borderaxespad=1)

    plt.show()

    if unmapped_frac_trimmed > unmapped_fraction_threshold:
        pass

    t5l, t3l, il = region_analysis(df, 3, 1)
    print("finished parsing")
    plt.hist(t5l, alpha=0.2, bins=np.arange(0,145,5), color='green')
    plt.xlim(0,145)
    #plt.axvline(x=61, linestyle='dashed', linewidth=2, color='red', alpha=0.8) # y-top of nanopore 1d
    #plt.axvline(x=120, linestyle='dashed', linewidth=2, color='red', alpha=0.8) # top of nanopore 1d2. 59 + 61
    plt.axvline(x=45, linestyle='dashed', linewidth=2, color='red', alpha=0.8) # pacbio
    plt.show()

    plt.hist(t3l, alpha=0.2, bins=np.arange(0, 145, 5), color='orange')
    plt.xlim(145,0)
    #plt.axvline(x=22, linestyle='dashed', linewidth=2, color='red', alpha=0.8) # y-bottom of nanopore 1d
    #plt.axvline(x=86, linestyle='dashed', linewidth=2, color='red', alpha=0.8) # bottom of nanopore 1d2. 22 + 64
    plt.axvline(x=45, linestyle='dashed', linewidth=2, color='red', alpha=0.8) # pacbio
    plt.show()

    """
    # internal break. experimental.
    plt.hist(il, alpha=0.2, color='blue')
    print(len(il))
    plt.plot(il, ([-0.5]*len(il)) + 0.1*np.random.random(len(il)), '+k')
    plt.show()
    """

    mid_threshold = 7 # ont
    mid_threshold = 8 # pb
    
    plt.grid(True)
    plt.boxplot([df[5].values[np.where(df[4] == 0.0)], df[5].values[np.where(df[4] != 0.0)]])
    plt.xticks([1,2], ["Non-OLP", "OLP"])
    ymin, ymax = plt.gca().get_ylim()
    plt.axhspan(0,  5, facecolor='red', alpha=0.1)
    plt.axhspan(5,  mid_threshold, facecolor='yellow', alpha=0.1)
    plt.axhspan(mid_threshold, ymax, facecolor='green', alpha=0.1)
    #plt.boxplot(df[5].values[np.where(df[4] == 0.0)])
    plt.ylim(0, ymax)
    plt.show()

    ### read score after size binning.
    subplot = gen_boxplot_length_vs_coverage(df, 3000.0)
    bin_size = df.groupby('Interval').size()
    boundary_reliable_bin = np.where(bin_size >= 50)[0].max()
    xmin, xmax = plt.gca().get_xlim()
    plt.axvspan(boundary_reliable_bin+1.5, xmax+1, facecolor='gray', alpha=0.1)
    plt.title("Read coverage over different length reads")
    plt.xticks(np.arange(xmax+1), [int(i) for i in np.arange(xmax+1)*3000])
    plt.ylim(0, mean_main+20*np.sqrt(cov_main))
    plt.suptitle("")
    #plt.savefig('Box_plot_quality.png')
    plt.show()
    plt.close()

    #interval = est_confidence_interval_params(df, 2, nbs=1000)
    #print(interval)
    #print(st.norm.ppf(0.2, loc=interval[0][1], scale=interval[1][0]))

