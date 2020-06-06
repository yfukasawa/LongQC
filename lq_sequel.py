import os, sys, logging, json, pysam
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot     as plt
import numpy                 as np
import xml.etree.ElementTree as et
import lq_gamma

from multiprocessing import Process
from operator        import itemgetter
from scipy.stats     import gaussian_kde, gamma
from lq_utils        import get_N50, get_NXX, rgb

# This module was inspired by and is translation of SEQUELstats provided by Vertebrate resequencing team
# github.com/VertebrateResequencing/SEQUELstats

def get_readtype(header):
    for d in header['RG']:
        if 'DS' in d and 'READTYPE' in d['DS']:
            vals = d['DS'].split(";")
            for v in vals:
                if v.split("=")[0] == 'READTYPE':
                    return v.split("=")[1]

def set_scrap(list, in_bam, snr):
    control_throughput = 0

    for r in in_bam:
        # scrap case
        if not r.has_tag('sz') or not r.has_tag('sc'):
            continue

        if r.get_tag('sz') == 'N':
            id_temp = r.query_name.split("/")
            zmw  = id_temp[1] # zmw
            _pos = id_temp[2].split("_")
            item  = (int(_pos[0]), int(_pos[1]), r.get_tag('sc')) # start, end
            if zmw not in list:
                list[zmw] = []
            list[zmw].append(item)

            #if r.has_tag('sn') and r.get_tag('sc') != "L":
            #    i = 0
            #    for f in r.get_tag('sn'):
            #        snr[i].append(f)
            #        i += 1

        elif r.get_tag('sz') == 'C':
            id_temp = r.query_name.split("/")
            zmw  = id_temp[1] # zmw
            _pos = id_temp[2].split("_")

            if r.get_tag('sc') == 'F':
                control_throughput += int(_pos[1]) - int(_pos[0]) + 1

    return control_throughput

def set_subreads(list, in_bam, snr):
    for r in in_bam:
        # scrap case
        id_temp = r.query_name.split("/")
        zmw  = id_temp[1] # zmw
        _pos = id_temp[2].split("_")
        pos  = (int(_pos[0]), int(_pos[1]), "S") # start, end           
        if zmw not in list:
            list[zmw] = []
        list[zmw].append(pos)

        if r.has_tag('sn'):
            i = 0
            for f in r.get_tag('sn'):
                snr[i].append(f)
                i += 1

# l: list of tuple, which contains start, end, and read class
def construct_polread(l):
    _end = 0
    _hs = _he =-1
    s_flag = a_flag = False

    tot = 0
    hq  = 0
    ad_num = 0

    ql_cigar_like = [] # quality
    st_cigar_like = [] # seq type
    
    s_l = sorted(l, key=itemgetter(0,1))

    for i in s_l:
        s = i[0]
        e = i[1]
        c = i[2]

        if _end != 0 and _end != s:
            if _hs >= 0:
                hq -= s-_end-1
            ql_cigar_like.append('%d%s' % (s-_end-1, "G"))
            st_cigar_like.append('%d%s' % (s-_end-1, "G"))
            tot += s-_end-1

        _end = e

        if c == 'L':
            if _hs >= 0:
                hq += _he - _hs
                ql_cigar_like.append('%d%s' % (_he-_hs+1, "H"))
                _he = _hs = -1
            ql_cigar_like.append('%d%s' % (e-s+1, c))
        else:
            if _hs < 0:
                _hs = s
            _he = e

            if c == 'S':
                s_flag = True
            elif c == 'A':
                a_flag = True
                ad_num += 1
            
        tot += e-s #old e and new s can be the same. do not add +1.
        st_cigar_like.append('%d%s' % (e-s+1, c))

    if _hs >= 0:
        hq += _he - _hs
        ql_cigar_like.append('%d%s' % (_he-_hs+1, "H"))

    if hq > 0:
        hq  += 1
    tot += 1
    
    #if s_flag and a_flag:
    if s_flag:
        # polymerase read
        return ("".join(ql_cigar_like), "".join(st_cigar_like), hq, tot, True, ad_num)
    else:
        return ("".join(ql_cigar_like), "".join(st_cigar_like), hq, tot, False, ad_num)

def parse_sts_xml(filepath, ns=None):
    tree = et.parse(filepath)
    root = tree.getroot()

    #print( root.find(".//{%s}HqRegionSnrDist[@Channel='A']" % ns) )

    bc = root.findall("./{http://pacificbiosciences.com/PacBioPipelineStats.xsd}ProdDist/{%s}BinCounts"  % ns)
    bl = root.findall("./{http://pacificbiosciences.com/PacBioPipelineStats.xsd}ProdDist/{%s}BinLabels"  % ns)
    
    p0 = p1 = p2 = 0
    for i,c in enumerate(bl[0]):
        if 'BinLabel' in c.tag:
            if 'Empty' in c.text:
                p0 = int(bc[0][i].text)
            elif 'Productive' in c.text:
                p1 = int(bc[0][i].text)
            elif 'Other' in c.text:
                p2 = int(bc[0][i].text)

    tot = p0 + p1 + p2
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

def get_bam_path(d, logger):
    subread_p = None
    scrap_p   = None

    if not os.path.isdir(d):
        logger.info("%s is not a dir" % d)
        return None

    list = os.listdir(d)

    for i in list:
        p = os.path.join(d, i)
        if os.path.isdir(p):
            pass
            continue
        if p.endswith(".scraps.bam"):
            scrap_p = p
            continue
        if p.endswith(".subreads.bam"):
            subread_p = p
            continue
                
    if subread_p and scrap_p:
        logger.info("Subreads bam file:%s" % subread_p)
        logger.info("Scraps bam file:%s" % scrap_p)
        return [subread_p, scrap_p]
    else:
        logger.ERROR("bam files are missing in %s" % d)
        return [None] * 2

def run_platformqc(data_path, output_path, *, suffix=None, b_width = 1000):
    if not suffix:
        suffix = ""
    else:
        suffix = "_" + suffix
    log_path     = os.path.join(output_path, "log", "log_sequel_platformqc" + suffix + ".txt")
    fig_path     = os.path.join(output_path, "fig", "fig_sequel_platformqc_length" + suffix + ".png")
    fig_path_bar = os.path.join(output_path, "fig", "fig_sequel_platformqc_adapter" + suffix + ".png")
    json_path    = os.path.join(output_path, "QC_vals_sequel"  + suffix + ".json")
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

    logger.info("Started sequel platform QC for %s" % data_path)

    # sequel
    xml_file = get_sts_xml_path(data_path, logger)

    if not xml_file:
        logger.warning("sts.xml is missing. Productivity won't be shown")
        [p0, p1, p2] = [None] * 3
    else:
        [p0, p1, p2] = parse_sts_xml(xml_file, ns="http://pacificbiosciences.com/PacBioBaseDataModel.xsd")
        logger.info("Parsed sts.xml")

    [subr_bam_p, scrap_bam_p] = get_bam_path(data_path, logger)
    if subr_bam_p and scrap_bam_p:
        scrap_bam = pysam.AlignmentFile(scrap_bam_p, 'rb', check_sq=False)
        subr_bam  = pysam.AlignmentFile(subr_bam_p, 'rb', check_sq=False)
    else:
        logger.ERROR("Platform QC failed due to missing bam files")
        return 1

    bam_reads          = {}
    snr                = [[], [], [], []]
    hr_fraction        = []
    tot_lengths        = []
    hr_lengths         = []
    ad_num_stat        = {}
    control_throughput = 0

    if get_readtype(scrap_bam.header) == 'SCRAP':
        logger.info("Started to load scraps.bam...")
        control_throughput = set_scrap(bam_reads, scrap_bam, snr)
    else:
        logger.ERROR("the given scrap file has incorrect header.")

    logger.info("Scrap reads were loaded.")

    if get_readtype(subr_bam.header) == 'SUBREAD':
        logger.info("Started to load subreads.bam...")
        set_subreads(bam_reads, subr_bam, snr)
    else:
        logger.ERROR("the given subread file has incorrect header.")

    logger.info("Subreads were loaded.")

    for k, v in bam_reads.items():
        #print(k)
        l = construct_polread(v)

        #print(l)
        if l[4]:
            hr_fraction.append(l[2]/l[3])
            tot_lengths.append(l[3])
            hr_lengths.append(l[2])
            if l[5] in ad_num_stat:
                ad_num_stat[l[5]] += 1
            else:
                ad_num_stat[l[5]] = 1


    max_adnum = max(ad_num_stat.keys())
    min_adnum = min(ad_num_stat.keys())

    left   = []
    height = [] 
    for i in range(min_adnum, max_adnum+1):
        left.append(i)
        if i in ad_num_stat:
            height.append(ad_num_stat[i])
        else:
            height.append(0)

    plt.bar(left, height)
    plt.savefig(fig_path_bar, bbox_inches="tight")
    plt.close()
    logger.info("Plotted bar plot for adpter occurence")

    (a, b) = lq_gamma.estimate_gamma_dist_scipy(hr_lengths)
    logger.info("Fitting by Gamma dist finished.")

    _max  = np.array(hr_lengths).max()
    _mean = np.array(hr_lengths).mean()
    _n50  = get_N50(hr_lengths)
    _n90  = get_NXX(hr_lengths, 90)
    throughput = np.sum(hr_lengths)
    longest    = np.max(hr_lengths)
    fracs      = np.mean(hr_fraction)

    tobe_json["Productivity"]         = {"P0": p0, "P1": p1, "P2":p2}
    tobe_json["Throughput"]           = int(throughput)
    tobe_json["Throughput(Control)"]  = int(control_throughput)
    tobe_json["Longest_read"]         = int(_max)
    tobe_json["Num_of_reads"]         = len(hr_lengths)
    tobe_json["polread_gamma_params"] = [float(a), float(b)]
    tobe_json["Mean_polread_length"]  = float(_mean)
    tobe_json["N50_polread_length"]   = float(_n50)
    tobe_json["Mean_HQ_fraction"]     = float(np.mean(fracs))
    tobe_json["Adapter_observation"]  = ad_num_stat

    with open(json_path, "w") as f:
        logger.info("Quality measurements were written into a JSON file: %s" % json_path)
        json.dump(tobe_json, f, indent=4)

    x = np.linspace(0, gamma.ppf(0.99, a, 0, b))
    est_dist = gamma(a, 0, b)
    plt.plot(x, est_dist.pdf(x), c=rgb(214,39,40) )
    plt.grid(True)
    plt.hist(hr_lengths, histtype='step', bins=np.arange(min(hr_lengths), _max + b_width, b_width), color=rgb(214,39,40), alpha=0.7, density=True)
    plt.xlabel('Read length')
    plt.ylabel('Probability density')

    if _mean >= 10000: # pol read mean is expected >= 10k and <= 15k, but omit the <= 15k condition.
        plt.axvline(x=_mean, linestyle='dashed', linewidth=2, color=rgb( 44, 160, 44), alpha=0.8)
    else:
        plt.axvline(x=_mean, linestyle='dashed', linewidth=2, color=rgb(188, 189, 34), alpha=0.8)

    if _n50 >= 20000:
        plt.axvline(x=_n50, linewidth=2, color=rgb( 44, 160, 44), alpha=0.8)
    else:
        plt.axvline(x=_n50, linewidth=2, color=rgb(188, 189, 34), alpha=0.8)

    plt.hist(tot_lengths, histtype='step', bins=np.arange(min(tot_lengths),max(tot_lengths) + b_width, b_width), color=rgb(31,119,180), alpha=0.7, density=True) 

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

    logger.info("Figs were generated.")
    logger.info("Finished all processes.")

# test
if __name__ == "__main__":
    # check_sq has to be False; otherwise, this doesn't work.
    run_platformqc("/home/fukasay/rawdata/pb/rs2_ecoli_pacbio_official/", "/home/fukasay/analyses/longQC/sequel_platform_test/")

    # fraction hist
    #plt.hist(hr_fraction, histtype='bar', bins=np.arange(0.0, 1.0 + 0.02, 0.02), color='blue') 
    #plt.show()

    """ # SNR PLOT
    fig, ax = plt.subplots(2,2, figsize=(6, 4))
    ax[0, 0].hist(snr[0], histtype='step', bins=np.arange(2,14 + 1, 1), color='blue' )
    ax[0, 1].hist(snr[1], histtype='step', bins=np.arange(2,14 + 1, 1), color='purple')
    ax[1, 0].hist(snr[2], histtype='step', bins=np.arange(2,14 + 1, 1), color='green')
    ax[1, 1].hist(snr[3], histtype='step', bins=np.arange(2,14 + 1, 1), color='red')

    # minimums are 5.5 and 4.0 for (A,C) and (G,T), respectively. But this is really true?
    #ls = np.linspace(2, 14)
    #kernel_A = gaussian_kde(df['SnrMean_A'].values[df['ReadScore']>0.1])
    #ax[0, 0].plot(ls, kernel_A(ls))
    #peak_A = ls[np.argmax(kernel_A(ls))]

    
    plt.show()
    """
