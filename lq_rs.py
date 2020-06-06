import pandas                as pd
import numpy                 as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot     as plt
import xml.etree.ElementTree as et
import os, logging, json
import lq_gamma
from scipy.stats import gaussian_kde, gamma
from lq_utils    import get_N50, get_NXX, rgb


def load_sts_csv(file_path):
    df = pd.read_table(file_path, sep=',')
    return df

def gen_boxplot_length_vs_score(df, interval):
    df['Interval'] = np.floor((df['HQRegionEnd'].values - df['HQRegionStart'].values)/interval)
    return df.boxplot(column='ReadScore', by='Interval', sym='+', rot=90, figsize=(int(max(df['Interval'])/5+0.5),6))

def gen_hist_tot_vs_hq_length(df):
    pass

def gen_snr_plots(df):
    # https://stackoverflow.com/questions/9767241/setting-a-relative-frequency-in-a-matplotlib-histogram
    fig, ax = plt.subplots(2,2, figsize=(6, 4))
    ax[0, 0].hist(df['SnrMean_A'].values[df['ReadScore']>0.1], histtype='step', bins=np.arange(2,14 + 0.1, 0.1), color='blue' )
    ax[0, 1].hist(df['SnrMean_C'].values[df['ReadScore']>0.1], histtype='step', bins=np.arange(2,14 + 0.1, 0.1), color='purple')
    ax[1, 0].hist(df['SnrMean_G'].values[df['ReadScore']>0.1], histtype='step', bins=np.arange(2,14 + 0.1, 0.1), color='green')
    ax[1, 1].hist(df['SnrMean_T'].values[df['ReadScore']>0.1], histtype='step', bins=np.arange(2,14 + 0.1, 0.1), color='red')

    # minimums are 5.5 and 4.0 for (A,C) and (G,T), respectively. But this is really true?
    ls = np.linspace(2, 14)
    kernel_A = gaussian_kde(df['SnrMean_A'].values[df['ReadScore']>0.1])
    ax[0, 0].plot(ls, kernel_A(ls))
    peak_A = ls[np.argmax(kernel_A(ls))]
  
    plt.show()

def parse_sts_xml(filepath, ns=None):
    tree = et.parse(filepath)
    root = tree.getroot()

    bc = root.findall("./{%s}ProdDist/{%s}BinCount" % (ns, ns))
    bl = root.findall("./{%s}ProdDist/{%s}BinLabel" % (ns, ns))

    p0 = p1 = p2 = 0

    for i,c in enumerate(bl):

        if 'BinLabel' in c.tag:
            if 'Empty' in c.text:
                p0 = int(bc[i].text)
            elif 'Productive' in c.text:
                p1 = int(bc[i].text)
            elif 'Other' in c.text:
                p2 = int(bc[i].text)

    return [p0, p1, p2]

def get_sts_xml_path(d, logger):
    if not os.path.isdir(d):
        logger.info("%s is not a dir" % d)
        return None

    list = os.listdir(d)

    for i in list:
        p = os.path.join(d, i)
        if os.path.isdir(p):
            pass
        if p.endswith(".sts.xml"):
            return p

    return None

def get_sts_csv_path(d, logger):
    if not os.path.isdir(d):
        logger.info("%s is not a dir" % d)
        return None

    list = os.listdir(d)

    for i in list:
        p = os.path.join(d, i)
        if os.path.isdir(p):
            pass
        if p.endswith(".sts.csv"):
            return p

    return None

def run_platformqc(data_path, output_path, *, suffix=None, b_width = 1000):
    if not suffix:
        suffix = ""
    else:
        suffix = "_" + suffix
    log_path  = os.path.join(output_path, "log", "log_rs2_platformqc" + suffix + ".txt")
    fig_path  = os.path.join(output_path, "fig", "fig_rs2_platformqc_length" + suffix + ".png")
    fig_path2 = os.path.join(output_path, "fig", "fig_rs2_platformqc_score" + suffix + ".png")
    json_path = os.path.join(output_path, "QC_vals_rs" + suffix + ".json")

    # json
    tobe_json = {}

    # output_path will be made too.
    if not os.path.isdir(os.path.join(output_path, "log")):
        os.makedirs(os.path.join(output_path, "log"), exist_ok=True)

    if not os.path.isdir(os.path.join(output_path, "fig")):
        os.makedirs(os.path.join(output_path, "fig"), exist_ok=True)

    ### logging conf ###
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    fh = logging.FileHandler(log_path, 'w')
    sh = logging.StreamHandler()

    formatter = logging.Formatter('%(module)s:%(asctime)s:%(lineno)d:%(levelname)s:%(message)s')
    fh.setFormatter(formatter)
    sh.setFormatter(formatter)

    logger.addHandler(sh)
    logger.addHandler(fh)
    #####################

    logger.info("Started RS-II platform QC for %s" % data_path)

    xml_file = get_sts_xml_path(data_path, logger)

    if not xml_file:
        logger.warning("sts.xml is missing. Productivity won't be shown")
        [p0, p1, p2] = [None] * 3
    else:
        [p0, p1, p2] = parse_sts_xml(xml_file, ns="http://pacificbiosciences.com/PipelineStats/PipeStats.xsd")
        logger.info("Parsed sts.xml")

    csv_path = get_sts_csv_path(data_path, logger)
    if not csv_path:
        logger.ERROR("Platform QC failed due to missing csv files")
        return 1

    df   = load_sts_csv(csv_path)
    logger.info("Stat file was loaded.")

    vals = df['HQRegionEnd'].values[df['ReadScore']>0.1] - df['HQRegionStart'].values[df['ReadScore']>0.1]

    (a, b) = lq_gamma.estimate_gamma_dist_scipy(vals, logger)
    logger.info("Fitting by Gamma dist finished.")
    _max   = np.array(vals).max()
    _mean  = np.array(vals).mean()
    _n50   = get_N50(vals)
    _n90   = get_NXX(vals, 90)
    throughput = np.sum(vals)

    ### HQ fraction over numbases
    fracs = vals/df['NumBases'].values[df['ReadScore']>0.1]
    #plt.hist(vals, histtype='bar', bins=np.arange(0.0, 1.0 + 0.02, 0.02), color='blue') 
    #plt.show()

    tobe_json["Productivity"]    = {"P0": p0, "P1": p1, "P2":p2}
    tobe_json["Throughput"]      = int(throughput)
    tobe_json["Longest_read"]    = int(_max)
    tobe_json["Num_of_reads"]    = len(vals)
    tobe_json["polread_gamma_params"] = [float(a), float(b)]
    tobe_json["Mean_polread_length"]  = float(_mean)
    tobe_json["N50_polread_length"]   = float(_n50)
    tobe_json["Mean_HQ_fraction"]     = float(np.mean(fracs))

    with open(json_path, "w") as f:
        logger.info("Quality measurements were written into a JSON file: %s" % json_path)
        json.dump(tobe_json, f, indent=4)

    x = np.linspace(0, gamma.ppf(0.99, a, 0, b))
    est_dist = gamma(a, 0, b)
    plt.plot(x, est_dist.pdf(x), c=rgb(214,39,40) )
    plt.grid(True)
    plt.hist(vals, histtype='step', bins=np.arange(min(vals), _max + b_width, b_width), color=rgb(214,39,40), alpha=0.7, density=True)
    plt.xlabel('Read length')
    plt.ylabel('Probability density')

    if _mean >= 10000: # pol read mean is expected >= 10k and <= 15k, but omit the <= 15k condition.
        plt.axvline(x=_mean, linestyle='dashed', linewidth=2, color=rgb( 44, 160, 44), alpha=0.8)
    else:
        plt.axvline(x=_mean, linestyle='dashed', linewidth=2, color=rgb(188, 189, 34), alpha=0.8)

    if _n50 >= 20000: # recent brochure says 20kb, but some old announcement says 14kb. let's see
        plt.axvline(x=_n50, linewidth=2, color=rgb( 44, 160, 44), alpha=0.8)
    else:
        plt.axvline(x=_n50, linewidth=2, color=rgb(188, 189, 34), alpha=0.8)

    vals = df['NumBases'].values[df['ReadScore']>0.1]
    plt.hist(vals, histtype='step', bins=np.arange(min(vals),max(vals) + b_width, b_width), color=rgb(31,119,180), alpha=0.7, density=True) 

    ymin, ymax = plt.gca().get_ylim()
    xmin, xmax = plt.gca().get_xlim()
    plt.text(xmax*0.6, ymax*0.72, r'$\alpha=%.3f,\ \beta=%.3f$' % (a,b) )
    plt.text(xmax*0.6, ymax*0.77,  r'Gamma dist params:' )

    plt.text(xmax*0.6, ymax*0.85, r'sample mean: %.3f' % (_mean,) )
    plt.text(xmax*0.6, ymax*0.9,  r'N50: %.3f' % (_n50,) )
    plt.text(xmax*0.6, ymax*0.95, r'N90: %.3f' % (_n90,) )

    plt.text(_mean, ymax*0.85, r'Mean')
    plt.text(_n50, ymax*0.9, r'N50')

    plt.savefig(fig_path, bbox_inches="tight")
    plt.close()
    #plt.show()

    ### read score after size binning.
    subplot = gen_boxplot_length_vs_score(df, b_width)
    xmin, xmax = plt.gca().get_xlim()
    plt.title("Read scores over different length reads")
    plt.xticks(np.arange(xmax+1), [int(i) for i in np.arange(xmax+1)*b_width])
    plt.suptitle("")
    plt.savefig(fig_path2, bbox_inches="tight")
    #plt.show()
    plt.close()

    logger.info("Figs were generated.")
    logger.info("Finished all processes.")


# test
if __name__ == "__main__":
    run_platformqc("/home/fukasay/basecalled/rs2/", "/home/fukasay/analyses/longQC/rs2_platform_test/")

    ### SNR
    #df = load_sts_csv('m170304_003258_42276_c101158722550000001823254607191760_s1_p0.sts.csv')
    #gen_snr_plots(df)
