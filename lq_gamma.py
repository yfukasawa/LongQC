from scipy import special
from scipy.stats import dgamma
from scipy.stats import gamma
from lq_utils    import rgb
from lq_utils    import rgb, get_N50
import sys
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from logging import getLogger
logger = getLogger(__name__)

def calc_ll(vals, a, b):
    return np.sum( np.log( dgamma.pdf(vals, a, loc=0, scale=1/b) ) )

# Estimating a Gamma distribution, T. Minka (2002)
# in case of many samples, (maybe) underflow occurred.
def estimate_gamma_dist_minka(vals, verbose=1):
    _v = np.array(vals)

    mean_v = np.mean(_v)

    a = 0.5 / ( np.log(mean_v) - np.mean( np.log(_v) ) )
    b = mean_v / a

    if verbose > 1:
        print("[Gamma Estimation] initial params are a = %f, b = %f." % (a,b))

    eps = float("inf")
    ll_p = calc_ll(_v, a, b)

    while eps > 0.000001:
        a = 1/a + ( np.mean( np.log( _v ) ) - special.polygamma(0, a) - np.log(mean_v) + np.log(a) )/( a**2 * (1/a - special.polygamma(1, a) ) )
        a = 1/a
        b = mean_v / a

        ll   = calc_ll(_v, a, b)
        eps  = np.abs( ll - ll_p )
        ll_p = ll

    return (a, b)

# simply call scipy method
def estimate_gamma_dist_scipy(vals):
    # shifting can happen (e.g. size selection of DNA fragments), but tentatively ignores this.
    alpha_hat, loc_hat, beta_hat = gamma.fit(vals, floc=0.0)

    logger.info("estimated Gamma dist params are a = %f, b = %f." % (alpha_hat, beta_hat))

    return (alpha_hat, beta_hat)


def plot_length_dist(fig_path, lengths, g_a, g_b, _max, _mean, _n50, isPb=False, b_width = 1000):
    x = np.linspace(0, gamma.ppf(0.99, g_a, 0, g_b))
    est_dist = gamma(g_a, 0, g_b)
    plt.hist(lengths, histtype='step', bins=np.arange(min(lengths),_max + b_width, b_width), color=rgb(214,39,40), alpha=0.7, density=True)
    plt.grid(True)
    plt.xlabel('Read length')
    plt.ylabel('Probability density')
    plt.axvline(x=_mean, linestyle='dashed', linewidth=2, color=rgb(214,39,40), alpha=0.8)
    plt.axvline(x=_n50,  linewidth=2, color=rgb(214,39,40), alpha=0.8)
    plt.xlim(0, gamma.ppf(0.99, g_a, 0, g_b))

    ymin, ymax = plt.gca().get_ylim()
    xmin, xmax = plt.gca().get_xlim()

    if not isPb:
        plt.text(xmax*0.6, ymax*0.72, r'$\alpha=%.3f,\ \beta=%.3f$' % (g_a, g_b) )
        plt.text(xmax*0.6, ymax*0.77,  r'Gamma dist params:' )
        plt.plot(x, est_dist.pdf(x), color=rgb(214,39,40) )

    plt.text(xmax*0.6, ymax*0.85, r'sample mean: %.3f' % (_mean,) )
    plt.text(xmax*0.6, ymax*0.9, r'N50: %.3f' % (_n50,) )

    plt.text(_mean, ymax*0.85, r'Mean', color=rgb(214,39,40))
    plt.text(_n50, ymax*0.9, r'N50', color=rgb(214,39,40))

    plt.axis('tight')
    plt.xlim(0, gamma.ppf(0.99, g_a, 0, g_b))
    plt.savefig(fig_path, bbox_inches="tight", transparent=True)
    #plt.show()
    plt.close()

# test
if __name__ == "__main__":
    fn  = sys.argv[1] # assumes output table of sdust in minimap2-coverage.
    out = sys.argv[2]

    df_mask = pd.read_table(fn, sep='\t', header=None)
    lengths = df_mask[2].values
    longest  = np.max(lengths)
    mean_len = np.array(lengths).mean()
    n50      = get_N50(lengths)

    (a, b) = estimate_gamma_dist_scipy(lengths)
    plt.rcParams['figure.figsize'] = (7, 7)
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    plot_length_dist(out, lengths, a, b, longest, mean_len, n50, False)
