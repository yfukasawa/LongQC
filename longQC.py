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
import pandas as pd

import lq_nanopore
import lq_rs
import lq_sequel

from time        import sleep
from scipy.stats import gamma
from lq_gamma    import estimate_gamma_dist_scipy
from lq_utils    import parse_fastq, parse_fasta, get_N50, rgb, sample_random_fastq, write_fastq
from lq_adapt    import cut_adapter
from lq_gcfrac   import plot_unmasked_gc_frac
from lq_exec     import LqExec
from lq_coverage import LqCoverage
from lq_mask     import LqMask

def command_run(args):
    if args.suf:
        suf = args.suf
    else:
        suf = None
    if args.platform == 'rs2':
        lq_rs.run_platformqc(args.raw_data_dir, args.out, suffix=suf)
    elif args.platform == 'sequel':
        lq_sequel.run_platformqc(args.raw_data_dir, args.out, suffix=suf)
    elif args.platform == 'minion':
        lq_nanopore.run_platformqc(args.platform, args.raw_data_dir, args.out, suffix=suf, n_channel=512)
    elif args.platform == 'gridion':
        lq_nanopore.run_platformqc(args.platform, args.raw_data_dir, args.out, suffix=suf, n_channel=512)
    else:
        pass

def command_help(args):
    print(parser.parse_args([args.command, '--help']))

def plot_length_dist(fig_path, lengths, g_a, g_b, _max, _mean, _n50, isPb=False, b_width = 1000):
    x = np.linspace(0, gamma.ppf(0.99, g_a, 0, g_b))
    est_dist = gamma(g_a, 0, g_b)
    plt.hist(lengths, histtype='step', bins=np.arange(min(lengths),_max + b_width, b_width), color=rgb(214,39,40), alpha=0.7, normed=True)
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
    plt.savefig(fig_path, bbox_inches="tight")
    #plt.show()
    plt.close()

def main(args):
    if hasattr(args, 'handler'):
        args.handler(args)
    else:
        parser.printhelp()

def command_sample(args):
    if args.suf:
        suffix = "_" + args.suf
    else:
        suffix = None

    cov_path    = os.path.join(args.out, "analysis", "minimap2", "coverage_out" + suffix + ".txt")
    cov_path_e  = os.path.join(args.out, "analysis", "minimap2", "coverage_err" + suffix + ".txt")
    sample_path = os.path.join(args.out, "analysis", "subsample" + suffix + ".fastq")
    log_path    = os.path.join(args.out, "log", "log_longQC_sampleqc" + suffix + ".txt")
    fig_path    = os.path.join(args.out, "fig", "fig_longQC_sampleqc_length" + suffix + ".png")
    fig_path_ma = os.path.join(args.out, "fig", "fig_longQC_sampleqc_masked_region" + suffix + ".png")
    fig_path_gc = os.path.join(args.out, "fig", "fig_longQC_sampleqc_gcfrac" + suffix + ".png")
    fig_path_cv = os.path.join(args.out, "fig", "fig_longQC_sampleqc_coverage" + suffix + ".png")
    fig_path_qv = os.path.join(args.out, "fig", "fig_longQC_sampleqc_olp_qv" + suffix + ".png")
    fig_path_ta = os.path.join(args.out, "fig", "fig_longQC_sampleqc_terminal_analysis" + suffix + ".png")
    fig_path_cl = os.path.join(args.out, "fig", "fig_longQC_sampleqc_coverage_over_read_length" + suffix + ".png")
    json_path   = os.path.join(args.out, "QC_vals_longQC_sampleqc" + suffix + ".json")

    df_mask = None

    # output_path will be made too.
    if not os.path.isdir(os.path.join(args.out, "analysis", "minimap2")):
        os.makedirs(os.path.join(args.out, "analysis", "minimap2"), exist_ok=True)

    if not os.path.isdir(os.path.join(args.out, "log")):
        os.makedirs(os.path.join(args.out, "log"), exist_ok=True)

    if not os.path.isdir(os.path.join(args.out, "fig")):
        os.makedirs(os.path.join(args.out, "fig"), exist_ok=True)

    ### logging conf ###
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
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

    logger.info("Computation of the low complexity region started.")
    lm = LqMask(reads, "/home/fukasay/Projects/minimap2_mod/sdust", args.out, suffix=suffix, n_proc=10)
    lm.run_async_sdust()
    logger.info("Summary table %s was made." % lm.get_outfile_path())
    lm.plot_masked_fraction(fig_path_ma)

    # list up seqs should be avoided
    df_mask      = pd.read_table(lm.get_outfile_path(), sep='\t', header=None)
    exclude_seqs = df_mask[(df_mask[2] > 500000) & (df_mask[3] > 0.2)][0].values.tolist() # len > 0.5M and mask_region > 20%. k = 15
    exclude_seqs = exclude_seqs + df_mask[(df_mask[2] > 20000) & (df_mask[3] > 0.4)][0].values.tolist() # len > 0.02M and mask_region > 40%. k = 12. more severe.
    logger.debug("Slowdown seqs list:\n%s" % "\n".join(exclude_seqs) )
    (sreads, s_n_seqs, s_n_bases) = sample_random_fastq(args.fastq, args.nsample, elist=exclude_seqs)
    logger.info('sequence sampling finished. #seqs:%d, #bases: %d, #n_sample: %f' % (s_n_seqs, s_n_bases, float(args.nsample)))
    write_fastq(sample_path, sreads)

    # asynchronized
    le = LqExec("/home/fukasay/Projects/minimap2_mod/minimap2-lite", logger=logger)
    le_args = shlex.split("%s -t %d -l 1000 %s %s" % (args.miniargs, int(args.minit), args.fastq, sample_path))
    le.exec(*le_args, out=cov_path, err=cov_path_e)

    logger.info("Overlap computation started. Process is %s" % le.get_pid())

    if df_mask is not None:
        lengths = df_mask[2].values
    else:
        lengths = []
    tobe_json = {}

    if args.adp5 and args.adp3:
        cut_adapter(reads, lengths, adp5, adp3, logger=logger)
    elif not args.adp5 and args.adp3:
        logger.warning("At present adapter processing using only 3\' is not supported.")
        #cut_adapter(reads, lengths, adp3, logger=logger)
    elif args.adp5 and not args.adp3:
        cut_adapter(reads, lengths, adp5, adp_b=None, logger=logger)

    if len(lengths) == 0:
        lengths = [len(r[1]) for r in reads]

    # length distribution. a ~= 1.0 is usual (exponential dist).
    (a, b) = estimate_gamma_dist_scipy(lengths, logger)
    
    throughput = np.sum(lengths)
    longest    = np.max(lengths)
    mean_len   = np.array(lengths).mean()
    n50        = get_N50(lengths)

    logger.info("Throughput: %d" % throughput)
    logger.info("Length of longest read: %d" % longest)
    logger.info("The number of reads: %d", len(lengths))

    tobe_json["Throughput"]       = int(throughput)
    tobe_json["Longest_read"]     = int(longest)
    tobe_json["Num_of_reads"]     = len(lengths)
    tobe_json["gamma_params"]     = [float(a), float(b)]
    tobe_json["Mean_read_length"] = float(mean_len)
    tobe_json["N50_read_length"]  = float(n50)

    plot_length_dist(fig_path, lengths, a, b, longest, mean_len, n50, True if args.pb else False)
    logger.info("Genarated the sample read length plot.")

    (gc_read_mean, gc_read_sd) = plot_unmasked_gc_frac(reads, logger=logger, fp=fig_path_gc)
    logger.info("Genarated the sample gc fraction plot.")

    tobe_json["Mean_GC_content"] = gc_read_mean
    tobe_json["SD_GC_content"]   = gc_read_sd

    # here wait until the minimap procerss finishes
    while True:
        if le.get_poll() is not None:
            logger.info("Process of %s terminated." % le.get_bin_path)
            break
        logger.info("Calculating overlaps of sampled reads...")
        sleep(30)

    logger.info("Overlap computation finished.")

    # execute minimap2_coverage
    logger.info("Generating coverage related plots...")
    lc = LqCoverage(cov_path, logger)
    lc.plot_coverage_dist(fig_path_cv)
    lc.plot_unmapped_frac_terminal(fig_path_ta)
    lc.plot_qscore_dist(fig_path_qv)
    lc.plot_length_vs_coverage(fig_path_cl)
    logger.info("Generated coverage related plots.")

    tobe_json["Non-overlapped fraction"] = float(lc.get_unmapped_frac())
    tobe_json["Mean_coverage"] = float(lc.get_mean())
    tobe_json["SD_coverage"]   = float(lc.get_sd())
    tobe_json["Estimated_crude_genome_size"] = str(lc.calc_genome_size(throughput))

    with open(json_path, "w") as f:
        logger.info("Quality measurements were written into a JSON file: %s" % json_path)
        json.dump(tobe_json, f, indent=4)

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
    platforms = ["rs2", "sequel", "minion", "gridion"]
    parser_run = subparsers.add_parser('runqc', help='see `runqc -h`')
    parser_run.add_argument('-s', '--suffix', help='suffix for each output file.', dest = 'suf', default = None)
    parser_run.add_argument('-o', '--output', help='path for output directory', dest = 'out', default = None)
    parser_run.add_argument('platform', choices=platforms, help='a platform to be evaluated. ['+", ".join(platforms)+']', metavar='platform')
    parser_run.add_argument('raw_data_dir', type=str, help='a path for a dir containing the raw data')
    #parser_run.add_argument('--rs', help='asseses a run of PacBio RS-II', dest = 'pbrs', action = 'store_true', default = None)
    #parser_run.add_argument('--sequel', help='asseses a run of PacBio Sequel', dest = 'pbsequel', choices=['kit2', 'kit2.1'], default = None)
    #parser_run.add_argument('--minion', help='asseses a run of ONT MinION', dest = 'ontmin', action = 'store_true', default = None)
    #parser_run.add_argument('--gridion', help='asseses a run of ONT GridION', dest = 'ontgrid', action = 'store_true', default = None)
    #parser_sample.add_argument('--promethion', help='asseses a run of ONT PromethION', dest = 'ontprom', action = 'store_true', default = None)
    parser_run.set_defaults(handler=command_run)

    # run sample
    adapter_sets = ["rs2", "sequel", "ont-ligation", "ont-rapid", "ont-1dsq"]
    parser_sample = subparsers.add_parser('sampleqc', help='see `sampleqc -h`')
    parser_sample.add_argument('--pb', help='sample data from PacBio sequencers', dest = 'pb', action = 'store_true', default = None)
    parser_sample.add_argument('--ont', help='sample data from ONT sequencers', dest = 'ont', action = 'store_true', default = None)
    parser_sample.add_argument('--adapter_set', help='adapter set. ['+", ".join(adapter_sets)+']', choices=adapter_sets, default = None, metavar='')
    parser_sample.add_argument('--adapter_5', help='adapter sequence for 5\'', dest = 'adp5', default = None)
    parser_sample.add_argument('--adapter_3', help='adapter sequence for 3\'', dest = 'adp3', default = None)
    parser_sample.add_argument('--n_sample', help='the number/fraction of sequences for sampling.', dest = 'nsample', default = 10000)
    parser_sample.add_argument('--minimap2_args', help='the arguments for minimap2.', dest = 'miniargs', default = '-Y -k 12 -w 5')
    parser_sample.add_argument('--minimap2_thread', help='the number of threads for sequences for minimap2.', dest = 'minit', default = 50)
    parser_sample.add_argument('-s', '--suffix', help='suffix for each output file.', dest = 'suf', default = None)
    parser_sample.add_argument('-o', '--output', help='path for output directory', dest = 'out', default = None)
    parser_sample.add_argument('fastq', help='Input in the Fastq format', type=str)
    parser_sample.set_defaults(handler=command_sample)

    # help
    parser_help = subparsers.add_parser('help', help='see `help -h`')
    parser_help.add_argument('command', help='')
    parser_help.set_defaults(handler=command_help)

    args = parser.parse_args()
    main(args)
