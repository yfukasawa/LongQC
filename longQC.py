'''
   Author: Yoshinori Fukasawa, Bioscience Core Lab @ KAUST, KSA

   Project Name: longQC.py
   Start Date: 2017-10-10
   Version: 0.1

   Usage:
      longQC.py [options]

      Try 'longQC.py -h' for more information.

    Purpose: longQC enables you to asses quality of sequence data
             coming from third-generation sequencers (long read).

    Bugs: Please contact yoshinori.fukasawa@kaust.edu.sa
'''

import sys, os, logging, json, argparse, shlex
import matplotlib.pyplot as plt
import numpy             as np
from scipy.stats import gamma
from lq_gamma    import estimate_gamma_dist_scipy
from lq_utils    import parse_fastq, parse_fasta, get_N50, rgb, sample_random_fastq, write_fastq
from lq_adapt    import cut_adapter
from lq_gcfrac   import plot_unmasked_gc_frac
from lq_exec     import LqExec
from lq_coverage import LqCoverage

def command_run(args):
    print(args)

def command_help(args):
    print(parser.parse_args([args.command, '--help']))

def plot_length_dist(fig_path, lengths, g_a, g_b, _max, _mean, _n50, b_width = 1000):
    
    x = np.linspace(0, gamma.ppf(0.99, g_a, 0, g_b))
    est_dist = gamma(g_a, 0, g_b)
    plt.hist(lengths, histtype='step', bins=np.arange(min(lengths),_max + b_width, b_width), color=rgb(214,39,40), alpha=0.7, normed=True)
    plt.plot(x, est_dist.pdf(x), color=rgb(214,39,40) )
    plt.grid(True)
    plt.xlabel('Read length')
    plt.ylabel('Probability density')
    plt.axvline(x=_mean, linestyle='dashed', linewidth=2, color=rgb(214,39,40), alpha=0.8)
    plt.axvline(x=_n50,  linewidth=2, color=rgb(214,39,40), alpha=0.8)
    plt.xlim(0, gamma.ppf(0.99, g_a, 0, g_b))

    ymin, ymax = plt.gca().get_ylim()
    xmin, xmax = plt.gca().get_xlim()

    plt.text(xmax*0.6, ymax*0.72, r'$\alpha=%.3f,\ \beta=%.3f$' % (g_a, g_b) )
    plt.text(xmax*0.6, ymax*0.77,  r'Gamma dist params:' )

    plt.text(xmax*0.6, ymax*0.85, r'sample mean: %.3f' % (_mean,) )
    plt.text(xmax*0.6, ymax*0.9, r'N50: %.3f' % (_n50,) )

    plt.text(_mean, ymax*0.85, r'Mean', color=rgb(214,39,40))
    plt.text(_n50, ymax*0.9, r'N50', color=rgb(214,39,40))

    plt.axis('tight')
    plt.xlim(0, gamma.ppf(0.99, g_a, 0, g_b))
    plt.savefig(fig_path, bbox_inches="tight")
    #plt.show()
    plt.close()

def main(args):
    if hasattr(args, 'handler'):
        args.handler(args)
    else:
        parser.print_help()

def command_sample(args):
    if args.suf:
        cov_path    = os.path.join(args.out, "analysis", "minimap2", "coverage_out_" + args.suf + ".txt")
        cov_path_e  = os.path.join(args.out, "analysis", "minimap2", "coverage_err_" + args.suf + ".txt")
        sample_path = os.path.join(args.out, "analysis", "subsample_" + args.suf + ".fastq")
        log_path    = os.path.join(args.out, "log", "log_longQC_sampleqc_" + args.suf + ".txt")
        fig_path    = os.path.join(args.out, "fig", "fig_longQC_sampleqc_length_" + args.suf + ".png")
        fig_path_gc = os.path.join(args.out, "fig", "fig_longQC_sampleqc_gcfrac_" + args.suf + ".png")
        json_path   = os.path.join(args.out, "QC_vals_longQC_sampleqc_" + args.suf + ".json")
    else:
        cov_path    = os.path.join(args.out, "analysis", "minimap2", "coverage_out.txt")
        cov_path_e  = os.path.join(args.out, "analysis", "minimap2", "coverage_err.txt")
        sample_path = os.path.join(args.out, "analysis", "subsample.fastq")
        log_path    = os.path.join(args.out, "log", "log_longQC_sampleqc.txt")
        fig_path    = os.path.join(args.out, "fig", "fig_longQC_sampleqc_length.png")
        fig_path_gc = os.path.join(args.out, "fig", "fig_longQC_sampleqc_gcfrac.png")
        json_path   = os.path.join(args.out, "QC_vals_longQC_sampleqc.json")

    # output_path will be made too.
    if not os.path.isdir(os.path.join(args.out, "analysis", "minimap2")):
        os.makedirs(os.path.join(args.out, "analysis", "minimap2"), exist_ok=True)

    if not os.path.isdir(os.path.join(args.out, "log")):
        os.makedirs(os.path.join(args.out, "log"), exist_ok=True)

    if not os.path.isdir(os.path.join(args.out, "fig")):
        os.makedirs(os.path.join(args.out, "fig"), exist_ok=True)

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

    (reads, n_seqs, n_bases) = parse_fastq(args.fastq)
    logger.info('fastq file parsing was finished. #seqs:%d, #bases: %d' % (n_seqs, n_bases))
    (sreads, s_n_seqs, s_n_bases) = sample_random_fastq(args.fastq, args.nsample)
    logger.info('sequence sampling finished. #seqs:%d, #bases: %d, #n_sample: %f' % (s_n_seqs, s_n_bases, float(args.nsample)))
    write_fastq(sample_path, sreads)

    # asynchronized
    le = LqExec("/home/fukasay/Projects/minimap2_mod/minimap2-lite", logger=logger)
    le_args = shlex.split("-Y -t 50 -l 0 -c 7 %s %s" % (sample_path, args.fastq))
    le.exec(*le_args, out=cov_path, err=cov_path_e)

    logger.info("Overlap computation started. Process is %d" % le.get_pid())

    lengths   = []
    tobe_json = {}

    cut_adapter(reads, lengths, logger, "AATGTACTTCGTTCAGTTACGTATTGCT")

    """
    print('traverse was finished. %d reads were processed.' % tot)
    print('Reads having adapter-like seq in both ends: %d / %d' % (trimboth, tot))
    print('Reads having adapter-like seq in 5\' ends: %d / %d' % (trim5 - trimboth, tot))
    print('Reads having adapter-like seq in 3\' ends: %d / %d' % (trim3 - trimboth, tot))
    print('Max identity in the reads: %.3f' % iden_max5)
    print('Max identity in the reads: %.3f' % iden_max3)
    print(max3_seq)
    """

    # length distribution. a ~= 1.0 is usual (exponential dist).
    (a, b) = estimate_gamma_dist_scipy(lengths, logger)
    
    throughput = np.sum(lengths)
    longest    = np.max(lengths)
    mean_len   = np.array(lengths).mean()
    n50        = get_N50(lengths)

    logger.info("Throughput: %d" % throughput)
    logger.info("Length of longest read: %d" % longest)
    logger.info("The number of reads: %d", len(lengths))

    #print("Throughput:", throughput, "Length of longest read:", longest, "The number of reads:", len(lengths))

    #(dip_statistic, pval) = diptest.diptest(lengths, boot_pval=True)
    #print (dip_statistic, pval)

    tobe_json["Throughput"]       = int(throughput)
    tobe_json["Longest_read"]     = int(longest)
    tobe_json["Num_of_reads"]     = len(lengths)
    tobe_json["gamma_params"]     = [float(a), float(b)]
    tobe_json["Mean_read_length"] = float(mean_len)
    tobe_json["N50_read_length"]  = float(n50)

    with open(json_path, "w") as f:
        logger.info("Quality measurements were written into a JSON file: %s" % json_path)
        json.dump(tobe_json, f, indent=4)

    plot_length_dist(fig_path, lengths, a, b, longest, mean_len, n50)
    logger.info("Genarated the sample read length plot.")

    plot_unmasked_gc_frac(fig_path_gc, reads, logger)
    logger.info("Genarated the sample gc fraction plot.")

    # here wait until the minimap procerss finishes
    while True:
        if le.get_poll() is not None:
            logger.info("Process of %s terminated." % le.get_bin_path)
            break

    logger.info("Overlap computation finished.")

    # execute minimap2_coverage
    lc = LqCoverage(cov_path, logger)
    lc.plot_coverage_dist()
    lc.plot_unmapped_frac_terminal()
    lc.plot_qscore_dist()
    lc.plot_length_vs_coverage()

    logger.info("Finished all processes.")

# stand alone
if __name__ == "__main__":
    # parsing
    parser = argparse.ArgumentParser(
        prog='LongQC.py',
        description='LongQC is a software to asses the quality of long read data from the third generation sequencers.',
        add_help=True,
    )
    subparsers = parser.add_subparsers()

    # run qc
    parser_run = subparsers.add_parser('runqc', help='see `runqc -h`')
    parser_run.add_argument('--rs', help='asseses a run of PacBio RS-II', dest = 'pbrs', action = 'store_true', default = None)
    parser_run.add_argument('--sequel', help='asseses a run of PacBio Sequel', dest = 'pbsequel', choices=['kit2', 'kit2.1'], default = None)
    parser_run.add_argument('--minion', help='asseses a run of ONT MinION', dest = 'ontmin', action = 'store_true', default = None)
    parser_run.add_argument('--gridion', help='asseses a run of ONT GridION', dest = 'ontgrid', action = 'store_true', default = None)
    #parser_sample.add_argument('--promethion', help='asseses a run of ONT PromethION', dest = 'ontprom', action = 'store_true', default = None)
    parser_run.set_defaults(handler=command_run)

    # run sample
    parser_sample = subparsers.add_parser('sampleqc', help='see `sampleqc -h`')
    parser_sample.add_argument('--pb', help='asseses a sample data from PacBio sequencers', dest = 'pb', action = 'store_true', default = None)
    parser_sample.add_argument('--ont', help='asseses a sample data from ONT sequencers', dest = 'ont', action = 'store_true', default = None)
    parser_sample.add_argument('--adapter_5', help='Specify adapter sequence for 5\'', dest = 'adp5', default = None)
    parser_sample.add_argument('--adapter_3', help='Specify adapter sequence for 3\'', dest = 'adp5', default = None)
    parser_sample.add_argument('--n_sample', help='Specify the number/fraction of sequences for sampling.', dest = 'nsample', default = 10000)
    parser_sample.add_argument('-s', '--suffix', help='Suffix for each output file.', dest = 'suf', default = None)
    parser_sample.add_argument('-o', '--output', help='Path for output directory', dest = 'out', default = None)
    parser_sample.add_argument('fastq', help='input in the Fastq format', type=str)
    parser_sample.set_defaults(handler=command_sample)

    # help
    parser_help = subparsers.add_parser('help', help='see `help -h`')
    parser_help.add_argument('command', help='')
    parser_help.set_defaults(handler=command_help)

    args = parser.parse_args()
    main(args)
