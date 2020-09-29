'''
   Author: Yoshinori Fukasawa, Bioscience Core Lab @ KAUST, KSA

   Project Name: longQC.py
   Start Date: 2017-10-10

   Usage:
      longQC.py [options]

      Try 'longQC.py -h' for more information.

    Purpose: LongQC enables you to asses the quality of sequence data
             coming from third-generation sequencers (long read).

    Bugs: Please contact to yoshinori.fukasawa@kaust.edu.sa
'''
import sys, os, json, argparse, shlex, array
import logging
import numpy  as np
import pandas as pd
from time        import sleep
from scipy.stats import gamma
from jinja2      import Environment, FileSystemLoader
from collections import OrderedDict
from multiprocessing import Pool
from _version import __version__

import lq_nanopore
import lq_rs
import lq_sequel

from lq_gamma    import estimate_gamma_dist_scipy, plot_length_dist
from lq_utils    import (eprint, open_seq_chunk, get_N50, subsample_from_chunk,
                         write_fastq, get_Qx_bases, copytree, guess_format, enc_b64_str)
from lq_adapt    import cut_adapter
from lq_gcfrac   import LqGC
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

def main(args):
    if hasattr(args, 'handler'):
        args.handler(args)
    else:
        parser.print_help()

def command_sample(args):
    if args.ncpu < 4:
        eprint("Error: -p/--ncpu needs to be 4 or higher.")
        sys.exit(1)

    if not os.path.exists(args.input):
        eprint("Error: input file %s does not exist." % args.input)
        sys.exit(1)

    if args.mem < 0 or args.mem > 2:
        eprint("Error: -m(--mem) option has an out-of-range ranged value.")
        sys.exit(1)

    if args.nsample < 0 or args.nsample > 10000:
        eprint("Error: -n(--n_sample) option has an out-of-range ranged value.")
        sys.exit(1)

    if os.path.exists(args.out):
        eprint("Error: output path %s already exists." % args.out)
        sys.exit(1)

    if args.db and args.short and args.ncpu < 9:
        print("-d/--db option with -b/--short flag requires at least nine cpus.")
        sys.exit(0)

    if args.db and args.ncpu < 6:
        print("-d/--db option requires at least six cpus.")
        sys.exit(0)

    if args.suf:
        suffix = "_" + args.suf
    else:
        suffix = ""

    ncpu = int(args.ncpu)
    path_minimap2  = os.path.join(os.path.dirname(os.path.abspath(__file__)), "minimap2-coverage")
    pb_control     = None
    merged_control = None
    cov_path    = os.path.join(args.out, "analysis", "minimap2", "coverage_out" + suffix + ".txt")
    cov_path_e  = os.path.join(args.out, "analysis", "minimap2", "coverage_err" + suffix + ".txt")
    sample_path = os.path.join(args.out, "analysis", "subsample" + suffix + ".fastq")
    if args.short:
        length_threshold = 500
        short_sample_path = os.path.join(args.out, "analysis", "short_subsample" + suffix + ".fastq")
        short_cov_path    = os.path.join(args.out, "analysis", "minimap2", "short_coverage_out" + suffix + ".txt")
        short_cov_path_e  = os.path.join(args.out, "analysis", "minimap2", "short_coverage_err" + suffix + ".txt")
        merged_cov_path   = os.path.join(args.out, "analysis", "minimap2", "merged_coverage_out" + suffix + ".txt")
    log_path    = os.path.join(args.out, "logs", "log_longQC_sampleqc" + suffix + ".txt")

    fig_path    = os.path.join(args.out, "figs", "fig_longQC_sampleqc_length" + suffix + ".png")
    fig_path_rq = os.path.join(args.out, "figs", "fig_longQC_sampleqc_average_qv" + suffix + ".png")
    fig_path_ma = os.path.join(args.out, "figs", "fig_longQC_sampleqc_masked_region" + suffix + ".png")
    fig_path_gc = os.path.join(args.out, "figs", "fig_longQC_sampleqc_gcfrac" + suffix + ".png")
    fig_path_cv = os.path.join(args.out, "figs", "fig_longQC_sampleqc_coverage" + suffix + ".png")
    fig_path_qv = os.path.join(args.out, "figs", "fig_longQC_sampleqc_olp_qv" + suffix + ".png")
    fig_path_ta = os.path.join(args.out, "figs", "fig_longQC_sampleqc_terminal_analysis" + suffix + ".png")
    fig_path_cl = os.path.join(args.out, "figs", "fig_longQC_sampleqc_coverage_over_read_length" + suffix + ".png")
    json_path   = os.path.join(args.out, "QC_vals_longQC_sampleqc" + suffix + ".json")
    fastx_path  = ""
    html_path   = os.path.join(args.out, "web_summary" + suffix + ".html")
    tempdb_path = ""

    df_mask = None
    minimap2_params = ''
    minimap2_db_params = ''
    minimap2_med_score_threshold = 0

    # for BCL
    '''
    db_index_pb  = os.path.join(os.path.dirname(os.path.abspath(__file__)), "db", "nt_pb.index")
    db_index_ont = os.path.join(os.path.dirname(os.path.abspath(__file__)), "db", "nt_ont.index")
    db_paf       = None
    db_paf_err   = None
    '''

    nonsense_read_error_threshold = 0.45
    nonsense_read_warn_threshold  = 0.25
    very_low_coverage_threshold   = 6

    # output_path will be made too.
    if not os.path.isdir(os.path.join(args.out, "analysis", "minimap2")):
        os.makedirs(os.path.join(args.out, "analysis", "minimap2"), exist_ok=True)

    if not os.path.isdir(os.path.join(args.out, "logs")):
        os.makedirs(os.path.join(args.out, "logs"), exist_ok=True)

    if not os.path.isdir(os.path.join(args.out, "figs")):
        os.makedirs(os.path.join(args.out, "figs"), exist_ok=True)

    ### logging conf ###
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    fh = logging.FileHandler(log_path, 'w')
    sh = logging.StreamHandler()

    formatter = logging.Formatter('%(module)s:%(asctime)s:%(lineno)d:%(levelname)s:%(message)s')
    fh.setFormatter(formatter)
    sh.setFormatter(formatter)

    logger.addHandler(sh)
    logger.addHandler(fh)
    #####################

    logger.info("Cmd: %s" % " ".join(sys.argv))

    if args.preset:
        p = args.preset
        if p == 'pb-rs2':
            args.pb = True
            args.adp5 = "ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT" if not args.adp5 else args.adp5
            args.adp3 = "ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT" if not args.adp3 else args.adp3
            minimap2_params    = "-Y -l 0 -q 160"
            minimap2_med_score_threshold = 80
            if args.short:
                minimap2_med_score_threshold_short = 60
        elif p == 'pb-sequel':
            args.pb = True
            args.sequel = True
            args.adp5 = "ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT" if not args.adp5 else args.adp5
            args.adp3 = "ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT" if not args.adp3 else args.adp3
            minimap2_params = "-Y -l 0 -q 160"
            minimap2_med_score_threshold = 80
            if args.short:
                minimap2_med_score_threshold_short = 60
        elif p == 'pb-hifi':
            args.pb = True
            args.sequel = True
            args.adp5 = "ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT" if not args.adp5 else args.adp5
            args.adp3 = "ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT" if not args.adp3 else args.adp3
            minimap2_params = "-Y -l 0 -q 160"
            minimap2_med_score_threshold = 80
        elif p == 'ont-ligation':
            args.ont = True
            args.adp5 = "AATGTACTTCGTTCAGTTACGTATTGCT" if not args.adp5 else args.adp5
            #args.adp3 = "GCAATACGTAACTGAACGAAGT"
            args.adp3 = "GCAATACGTAACTGAACG" if not args.adp3 else args.adp3
            minimap2_params = "-Y -l 0 -q 160"
            minimap2_med_score_threshold = 160
            if args.short:
                minimap2_med_score_threshold_short = 140
        elif p == 'ont-rapid':
            args.ont = True
            args.adp5 = "GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA" if not args.adp5 else args.adp5
            minimap2_params = "-Y -l 0 -q 160"
            minimap2_med_score_threshold = 160
            if args.short:
                minimap2_med_score_threshold_short = 140
        elif p == 'ont-1dsq':
            args.ont = True
            args.adp5 = "GGCGTCTGCTTGGGTGTTTAACCTTTTTGTCAGAGAGGTTCCAAGTCAGAGAGGTTCCT" if not args.adp5 else args.adp5
            args.adp3 = "GGAACCTCTCTGACTTGGAACCTCTCTGACAAAAAGGTTAAACACCCAAGCAGACGCCAGCAAT" if not args.adp3 else args.adp3
            minimap2_params = "-Y -l 0 -q 160"
            minimap2_med_score_threshold = 160
            if args.short:
                minimap2_med_score_threshold_short = 140

        if p == 'pb-hifi':
            if args.fast:
                minimap2_db_params = "-k 19 -w 10 -I %s" % args.inds
            else:
                minimap2_db_params = "-k 15 -w 5 -I %s" % args.inds
        else:
            if args.fast:
                minimap2_db_params = "-k 15 -w 5 -I %s" % args.inds
            else:
                minimap2_db_params = "-k 12 -w 5 -I %s" % args.inds 

        logger.info("Preset \"%s\" was applied. Options --pb(--ont) is overwritten." % (p,))

    file_format_code = guess_format(args.input)
    if file_format_code == 0:
        fastx_path = os.path.join(args.out, "analysis", "pbbam_converted_seq_file" + suffix + ".fastq")
        logger.info('Temporary work file was made at %s' % fastx_path)
    elif file_format_code == 4: #fast5
        fastx_path = os.path.join(args.out, "analysis", "fast5_converted_seq_file" + suffix + ".fastq")
        logger.info('Temporary work file was made at %s' % fastx_path)
    elif file_format_code == -1 or file_format_code == 1:
        logger.error('Input file is unsupported file format: %s' % args.input)
        sys.exit()
    else:
        fastx_path = args.input

    if args.pb:
        if args.sequel:
            filter_ref     = os.path.join(os.path.join(os.path.dirname(os.path.abspath(__file__)), "refs"), "Sequel_control_reference.fasta")
        else:
            filter_ref     = os.path.join(os.path.join(os.path.dirname(os.path.abspath(__file__)), "refs"), "RS2_control_reference.fasta")
        pb_control     = os.path.join(args.out, "analysis", "minimap2", "spiked_in_control" + suffix + ".txt")
        pb_control_err = os.path.join(args.out, "analysis", "minimap2", "spiked_in_control" + suffix + "_stderr.txt")
        minimap2_filtering_params = "-Y -Hk15 -w 10 -c 1 -l 0 --filter"
        nonsense_read_error_threshold = 0.2
        nonsense_read_warn_threshold  = 0.15
        if args.short:
            pb_control_short     = os.path.join(args.out, "analysis", "minimap2", "short_spiked_in_control" + suffix + ".txt")
            pb_control_short_err = os.path.join(args.out, "analysis", "minimap2", "short_spiked_in_control" + suffix + "_stderr.txt")
            merged_control       = os.path.join(args.out, "analysis", "minimap2", "merged_spiked_in_control" + suffix + ".txt")

    if args.short:
        minimap2_db_params_short = "-k 12 -w 5 -I %s" % args.inds 

    if args.db and file_format_code != 0:
        ncpu -= 3 # subtract cpus for the minimap2 db
        tempdb_path = os.path.join(args.out, "analysis", "minimap2", "t_db_minimap2")
        le = LqExec(os.path.join(path_minimap2, "minimap2-coverage"))
        le_args = shlex.split("%s -d %s %s" % (minimap2_db_params, tempdb_path, fastx_path))
        le.exec(*le_args, out=cov_path, err=cov_path_e)
        if args.short:
            ncpu -= 3 # subtract cpus further for the minimap2 db
            tempdb_short_path = os.path.join(args.out, "analysis", "minimap2", "t_db_minimap2_short")
            le_short = LqExec(os.path.join(path_minimap2, "minimap2-coverage"))
            le_args_short = shlex.split("%s -d %s %s" % (minimap2_db_params_short, tempdb_short_path, fastx_path))
            le_short.exec(*le_args_short, out=short_cov_path, err=short_cov_path_e)

    ### initialization for chunked reads ###
    pool = Pool(processes=2)
    ncpu -= 2 # subtract cpus for the executor
    pool_res  = {}
    lm = LqMask(os.path.join(path_minimap2, "sdust"), args.out, suffix=suffix, max_n_proc=10 if ncpu > 10 else ncpu)
    lg = LqGC(chunk_size=150)
    if args.adp5:
        num_trim5     = 0
        max_iden_adp5 = 0.0
        adp_pos5      = array.array('i')
    if args.adp3:
        num_trim3     = 0
        max_iden_adp3 = 0.0
        adp_pos3      = array.array('i')

    # vars for subsampling
    cum_n_seq = 0
    s_reads   = []
    #sample_random_fastq_list(reads, args.nsample, elist=exclude_seqs)
    chunk_n = 0
    for (reads, n_seqs, n_bases) in open_seq_chunk(args.input, file_format_code, chunk_size=args.mem*1024**3, is_upper=True):
        ### iterate over chunks
        ### 1. bam/fast5 to fastq conversion
        if file_format_code == 0 or file_format_code == 4:
            write_fastq(fastx_path, reads, is_chunk=True)

        ### 2. low-complexity region calc -> another process
        logger.info("Computation of the low complexity region started for a chunk %d" % chunk_n)
        lm.submit_sdust(reads, chunk_n)

        ### 3. adapter search -> another process
        if args.adp5 or args.adp3:
            logger.info("Adapter search is starting for a chunk %d." % chunk_n)
        if args.adp5 and args.adp3:
            #(tuple_5, tuple_3) = cut_adapter(reads, adp_t=args.adp5, adp_b=args.adp3, logger=logger)
            pool_res['adapter'] = pool.apply_async(cut_adapter, args=(reads,), kwds={'adp_t':args.adp5, 'adp_b':args.adp3})
        elif not args.adp5 and args.adp3:
            #tuple_3 = cut_adapter(reads, adp_b=args.adp3, adp_t=None, logger=logger)
            pool_res['adapter'] =  pool.apply_async(cut_adapter, args=(reads,), kwds={'adp_b':args.adp3})
        elif args.adp5 and not args.adp3:
            #tuple_5 = cut_adapter(reads, adp_t=args.adp5, adp_b=None, logger=logger)
            pool_res['adapter'] = pool.apply_async(cut_adapter, args=(reads,), kwds={'adp_t':args.adp5})

        ### 4. subsampling -> another process
        pool_res['subsample'] = pool.apply_async(subsample_from_chunk, args=(reads, cum_n_seq, s_reads, args.nsample))
        #pool_res['subsample'] = executor.submit(subsample_from_chunk, reads, cum_n_seq, s_reads, args.nsample, **{'minlen': 300})

        ### 5. GC fraction -> within this process as this is not pickable (class method)
        logger.info("Computation of the GC fraction started for a chunk %d" % chunk_n)
        lg.calc_read_and_chunk_gc_frac(reads)

        if args.adp5 and args.adp3:
            (tuple_5, tuple_3) = pool_res['adapter'].get()
            logger.info("Adapter search has done for a chunk %d." % chunk_n)
        elif not args.adp5 and args.adp3:
            tuple_3 = pool_res['adapter'].get()
            logger.info("Adapter search has done for a chunk %d." % chunk_n)
        elif args.adp5 and not args.adp3:
            tuple_5 = pool_res['adapter'].get()
            logger.info("Adapter search has done for a chunk %d." % chunk_n)

        ### 6. termination of one chunk
        s_reads = pool_res['subsample'].get()
        logger.info('subsample finished for chunk %d.' % chunk_n)

        # trimmed reads by edlib are saved as fastq
        if args.trim:
            write_fastq(args.trim, reads, is_chunk=True)
            logger.info("Trimmed read added.")
        if args.adp5 and tuple_5:
            if tuple_5[0] > max_iden_adp5:
                max_iden_adp5 = tuple_5[0]
            num_trim5 += tuple_5[1]
            adp_pos5.fromlist(tuple_5[2])
        if args.adp3 and tuple_3:
            if tuple_3[0] > max_iden_adp3:
                max_iden_adp3 = tuple_3[0]
            num_trim3 += tuple_3[1]
            adp_pos3.fromlist(tuple_3[2])

        chunk_n += 1
        cum_n_seq += n_seqs
    ### file traverse is over now.
    logger.info('Input file parsing was finished. #seqs:%d, #bases: %d' % (n_seqs, n_bases))

    # wait for completion of DUST analysis
    lm.close_pool()
    logger.info("Summary table %s was made." % lm.get_outfile_path())

    # list up seqs should be avoided
    df_mask      = pd.read_table(lm.get_outfile_path(), sep='\t', header=None)
    exclude_seqs = df_mask[(df_mask[2] > 500000) & (df_mask[3] > 0.2)][0].values.tolist() # len > 0.5M and mask_region > 20%. k = 15
    exclude_seqs = exclude_seqs + df_mask[(df_mask[2] > 10000) & (df_mask[3] > 0.4)][0].values.tolist() # len > 0.01M and mask_region > 40%. k = 12. more severe.
    logger.debug("Highly masked seq list:\n%s" % "\n".join(exclude_seqs) )

    # polishing subsampled seqs
    s_reads = [i for i in s_reads if i != 0] # removing empty pos. this happens if numseq < numsample
    ng_set  = set(exclude_seqs)
    ng_ovlp = 0
    ng_ovlp_indices = []
    for i, r in enumerate(s_reads):
        if r[0] in ng_set:
            ng_ovlp += 1
            ng_ovlp_indices.append(i)

    if ng_ovlp > 0:
        logger.info('There are %d overlap reads between highly masked samples and subsampled reads. Start replacing.' % ng_ovlp)
        temp = [0] * ng_ovlp
        j = 0
        for r in s_reads:
            ng_set.add(r[0]) # as skip already picked up ones
        for (reads, n_seqs, n_bases) in open_seq_chunk(args.input, file_format_code, chunk_size=0.1*1024**3):
            subsample_from_chunk(reads, j, temp, ng_ovlp, elist=ng_set)
            j += n_seqs
            if len([i for i in temp if i]) < ng_ovlp:
                continue
            else:
                break
        if len([i for i in temp if i]) < ng_ovlp:
            # an edgy case, but can happen.
            logger.warn('Replacing failed. Just removing highly masked ones.')
            for i in ng_ovlp_indices:
                s_reads[i] = 0
            s_reads = [i for i in s_reads if i]
        else:
            for i, t in enumerate(temp):
                logger.info('Replacing %s with %s.' % (s_reads[ng_ovlp_indices[i]][0], t[0]))
                s_reads[ng_ovlp_indices[i]] = t # replacing bad ones with ok ones

    s_n_seqs = len([i for i in s_reads if i])
    if args.short:
        ss_reads = [s for s in s_reads if len(s[1]) < length_threshold]
        if write_fastq(short_sample_path, ss_reads):
            logger.info('Short subsampled seqs were written to a file. #seqs:%d' % s_n_seqs)
        s_reads  = [s for s in s_reads if len(s[1]) >= length_threshold]
        if write_fastq(sample_path, s_reads):
            logger.info('Subsampled seqs were written to a file. #seqs:%d' % s_n_seqs)
    else:
        if write_fastq(sample_path, s_reads):
            logger.info('Subsampled seqs were written to a file. #seqs:%d' % s_n_seqs)

    # waiting db make by minimap2
    if args.db and file_format_code != 0:
        while True:
            if le.get_poll() is not None:
                logger.info("Process %s for %s terminated." % (le.get_pid(), le.get_bin_path()))
                break
            logger.info("Making a db of sampled reads...")
            sleep(10)
        if args.short:
            while True:
                if le_short.get_poll() is not None:
                    logger.info("Process %s for %s terminated." % (le_short.get_pid(), le_short.get_bin_path()))
                    break
                logger.info("Making a db of sampled short reads...")
                sleep(10)
        logger.info("Temp db %s was generated." % tempdb_path)

    # asynchronized minimap2 starts
    le = LqExec(os.path.join(path_minimap2, "minimap2-coverage"))
    if args.db and file_format_code != 0:
        le_args = shlex.split("%s -p %d -t %d %s %s" \
                              % (minimap2_params, int(minimap2_med_score_threshold), int(args.ncpu), tempdb_path, sample_path))
    else:
        le_args = shlex.split("%s %s -p %d -t %d %s %s" \
                              % (minimap2_params, minimap2_db_params, int(minimap2_med_score_threshold), int(args.ncpu), fastx_path, sample_path))
    le.exec(*le_args, out=cov_path, err=cov_path_e)
    logger.info("Overlap computation started. Process is %s" % le.get_pid())

    # gc frac plot
    gc_read_mean, gc_read_sd = lg.plot_unmasked_gc_frac(fp=fig_path_gc)
    logger.info("Genarated the sample gc fraction plot.")

    q7 = np.sum(df_mask[5].values) # make c code to compute Q7 now for speed
    #q10 =  get_Qx_bases(reads, threshold=10) # too slow
    logger.info("Q%d bases %d" % (7, q7))

    if df_mask is not None:
        lengths = df_mask[2].values
    else:
        logger.error("The reads summary table made by sdust does not exist!")
        sys.exit(1)

    tobe_json = {}
    
    # reads does not exist anymore due to chunking
    #if len(lengths) == 0:
    #    lengths = [len(r[1]) for r in reads]
    
    throughput = np.sum(lengths)
    longest    = np.max(lengths)
    mean_len   = np.array(lengths).mean()
    n50        = get_N50(lengths)

    # exceptionally short case.
    #if args.ont:
    #    if n50 < 1000 or float(len(np.where(np.asarray(lengths)< 1000)[0]))/len(lengths) > 0.25:
    #        minimap2_med_score_threshold = 60

    if n50 < 3000:
        lm.plot_qscore_dist(df_mask, 4, 2, interval=n50/2, fp=fig_path_rq)
    else:
        lm.plot_qscore_dist(df_mask, 4, 2, fp=fig_path_rq)

    # plot masked fraction
    lm.plot_masked_fraction(fig_path_ma)

    # length distribution. a ~= 1.0 is usual (exponential dist).
    (a, b) = estimate_gamma_dist_scipy(lengths)
    plot_length_dist(fig_path, lengths, a, b, longest, mean_len, n50, True if args.pb else False)
    logger.info("Genarated the sample read length plot.")

    logger.info("Throughput: %d" % throughput)
    logger.info("Length of longest read: %d" % longest)
    logger.info("The number of reads: %d", len(lengths))

    tobe_json["Yield"]            = int(throughput)
    tobe_json["Q7 bases"]         = str("%.2f%%" % float(100*q7/throughput))
    tobe_json["Longest_read"]     = int(longest)
    tobe_json["Num_of_reads"]     = len(lengths)
    tobe_json["Length_stats"] = {}
    tobe_json["Length_stats"]["gamma_params"]     = [float(a), float(b)]
    tobe_json["Length_stats"]["Mean_read_length"] = float(mean_len)
    tobe_json["Length_stats"]["N50_read_length"]  = float(n50)

    tobe_json["GC_stats"] = {}
    tobe_json["GC_stats"]["Mean_GC_content"] = float(gc_read_mean)
    tobe_json["GC_stats"]["SD_GC_content"]   = float(gc_read_sd)

    if args.adp5 and max_iden_adp5 >= 0.75:
        tobe_json["Stats_for_adapter5"] = {}
        tobe_json["Stats_for_adapter5"]["Num_of_trimmed_reads_5"] = num_trim5
        tobe_json["Stats_for_adapter5"]["Max_identity_adp5"] = max_iden_adp5
        tobe_json["Stats_for_adapter5"]["Average_position_from_5_end"] = np.mean(adp_pos5)
    if args.adp3 and max_iden_adp3 >= 0.75:
        tobe_json["Stats_for_adapter3"] = {}
        tobe_json["Stats_for_adapter3"]["Num_of_trimmed_reads_3"] = num_trim3
        tobe_json["Stats_for_adapter3"]["Max_identity_adp3"] = max_iden_adp3
        tobe_json["Stats_for_adapter3"]["Average_position_from_3_end"] = np.mean(adp_pos3)

    # here wait until the minimap procerss finishes
    while True:
        if le.get_poll() is not None:
            logger.info("Process %s for %s terminated." % (le.get_pid(), le.get_bin_path()))
            break
        logger.info("Calculating overlaps of sampled reads...")
        sleep(10)
    logger.info("Overlap computation finished.")

    if args.short:
        le_short = LqExec(os.path.join(path_minimap2, "minimap2-coverage"))
        if args.db and file_format_code != 0:
            le_short_args = shlex.split("%s -p %d -t %d %s %s" \
                                  % (minimap2_params, int(minimap2_med_score_threshold_short), int(args.ncpu), tempdb_short_path, short_sample_path))
        else:
            le_short_args = shlex.split("%s %s -p %d -t %d %s %s" \
                                  % (minimap2_params, minimap2_db_params_short, int(minimap2_med_score_threshold_short), int(args.ncpu), fastx_path, short_sample_path))
        le_short.exec(*le_short_args, out=short_cov_path, err=short_cov_path_e)
        logger.info("Overlap computation for short reads started. Process is %s" % le.get_pid())
        while True:
            if le_short.get_poll() is not None:
                logger.info("Process %s for %s terminated." % (le.get_pid(), le.get_bin_path()))
                break
            logger.info("Calculating overlaps of short sampled reads...")
            sleep(10)
        logger.info("Overlap computation for short reads finished.")
        with open(merged_cov_path, 'w') as outf:
            with open(cov_path, 'r') as inf:
                outf.write(inf.read())
            with open(short_cov_path, 'r') as inf:
                outf.write(inf.read())
        logger.info("Outputs for normal and short reads were merged.")
        
    # filtering for spiked in
    if args.pb:
        le_spike = LqExec(os.path.join(path_minimap2, "minimap2-coverage"))
        le_spike_args = shlex.split("%s -t %d %s %s" \
                                    % (minimap2_filtering_params, int(args.ncpu), filter_ref, sample_path))
        le_spike.exec(*le_spike_args, out=pb_control, err=pb_control_err)
        logger.info("Spike-in control filteration started. Process is %s" % le_spike.get_pid())

        # here wait until the minimap procerss finishes
        while True:
            if le.get_poll() is not None:
                logger.info("Process %s for %s terminated." % (le_spike.get_pid(), le_spike.get_bin_path()))
                break
            logger.info("Filtering spike-in control in sampled reads...")
            sleep(10)

        if args.short:
            le_spike_short = LqExec(os.path.join(path_minimap2, "minimap2-coverage"))
            le_spike_short_args = shlex.split("%s -t %d %s %s" \
                                              % (minimap2_filtering_params, int(args.ncpu), filter_ref, short_sample_path))
            le_spike_short.exec(*le_spike_short_args, out=pb_control_short, err=pb_control_short_err)
            logger.info("Spike-in control filteration started. Process is %s" % le_spike_short.get_pid())

            # here wait until the minimap procerss finishes
            while True:
                if le.get_poll() is not None:
                    logger.info("Process %s for %s terminated." % (le_spike_short.get_pid(), le_spike_short.get_bin_path()))
                    break
                logger.info("Filtering spike-in control in sampled reads...")
                sleep(10)
            logger.info("Filteration finished.")

            sleep(10)

            with open(merged_control, 'w') as outf:
                with open(pb_control, 'r') as inf:
                    outf.write(inf.read())
                with open(pb_control_short, 'r') as inf:
                    outf.write(inf.read())
            logger.info("Outputs for normal and short reads were merged.")
        logger.info("Filteration finished.")

    # for laggy file system, we neeed to wait a bit. otherwise, no data exception will be raised.
    sleep(10)

    # execute minimap2_coverage
    logger.info("Generating coverage related plots...")
    if args.short:
        lc = LqCoverage(merged_cov_path, isTranscript=args.transcript, control_filtering=merged_control)
    else:
        lc = LqCoverage(cov_path, isTranscript=args.transcript, control_filtering=pb_control)
    lc.plot_coverage_dist(fig_path_cv)
    lc.plot_unmapped_frac_terminal(fig_path_ta, \
                                   adp5_pos=np.mean(adp_pos5) if args.adp5 and adp_pos5 and np.mean(adp_pos5) > 0 else None, \
                                   adp3_pos=np.mean(adp_pos3) if args.adp3 and adp_pos3 and np.mean(adp_pos3) > 0 else None)
    lc.plot_qscore_dist(fig_path_qv)
    if n50 < 3000:
        lc.plot_length_vs_coverage(fig_path_cl, interval=n50/2)
    else:
        lc.plot_length_vs_coverage(fig_path_cl)
    logger.info("Generated coverage related plots.")

    if lc.is_no_coverage():
        pass # please keep this to avoid evaluate below conditions
    elif (args.transcript and float(lc.get_logn_mode()) < very_low_coverage_threshold) \
       or (lc.is_low_coverage() and float(lc.get_logn_mode()) < very_low_coverage_threshold) \
       or (float(lc.get_mean()) < very_low_coverage_threshold):
        logger.info("Coverage looks to be very low. Turns on the very low coverage mode.")
        very_low_coverage_mode = True

        if args.pb:
            nonsense_read_error_threshold = 0.1
            nonsense_read_warn_threshold  = 0.075
            
        '''
        # for BCL use. very low coverage mode
        le_db      = LqExec(os.path.join(path_minimap2, "minimap2"))
        db_paf     = os.path.join(args.out, "analysis", "minimap2", "db_hits" + suffix + ".paf")
        db_paf_err = os.path.join(args.out, "analysis", "minimap2", "db_hits_err" + suffix + ".txt")

        if args.pb:
            le_db_args = shlex.split("%s -t %d %s %s" \
                                     % ("-x map-pb --secondary=no", int(args.ncpu), db_index_pb, sample_path))
        elif args.ont:
            le_db_args = shlex.split("%s -t %d %s %s" \
                                     % ("-x map-ont --secondary=no", int(args.ncpu), db_index_ont, sample_path))
        le_db.exec(*le_db_args, out=db_paf, err=db_paf_err)
        logger.info("Very low coverage mode: DB retrieval has started. Process is %s" % le_db.get_pid())

        # here wait until the minimap procerss finishes
        while True:
            if le_db.get_poll() is not None:
                logger.info("Process %s for %s terminated." % (le_db.get_pid(), le_db.get_bin_path()))
                break
            logger.info("Very low coverage mode: DB retrieval is still on going...")
            sleep(100)

        lp = LqPaf(db_paf)
        mapped_ids = lp.get_mapped_reads_ids()
        '''
    else:
        # enough coverage to estimate stats
        very_low_coverage_mode = False

    # for laggy file system, we neeed to wait a bit. otherwise, no data exception will be raised.
    sleep(10)

    tobe_json["Coverage_stats"] = {}
    tobe_json["Coverage_stats"]["Estimated non-sense read fraction"] = float(lc.get_unmapped_med_frac())
    #tobe_json["Coverage_stats"]["Reliable Highly diverged fraction"] = float(lc.get_high_div_frac())
    if lc.get_control_frac():
        tobe_json["Coverage_stats"]['Estimated spiked-in control read fraction'] = float(lc.get_control_frac())

    if args.transcript:
        tobe_json["Coverage_stats"]["Mode_coverage"]  = float(lc.get_logn_mode())
        tobe_json["Coverage_stats"]["mu_coverage"]    = float(lc.get_logn_mu())
        tobe_json["Coverage_stats"]["sigma_coverage"] = float(lc.get_logn_sigma())
    elif lc.is_low_coverage():
        tobe_json["Coverage_stats"]["Mode_coverage"]  = float(lc.get_logn_mode())
        tobe_json["Coverage_stats"]["mu_coverage"]    = float(lc.get_logn_mu())
        tobe_json["Coverage_stats"]["sigma_coverage"] = float(lc.get_logn_sigma())
    elif lc.is_no_coverage():
        tobe_json["Coverage_stats"]["Mean_coverage"] = "NA"
        tobe_json["Coverage_stats"]["SD_coverage"]   = "NA"
    else:
        tobe_json["Coverage_stats"]["Mean_coverage"] = float(lc.get_mean())
        tobe_json["Coverage_stats"]["SD_coverage"]   = float(lc.get_sd())
        # adjust threshold for very low coverage

    tobe_json["Coverage_stats"]["Estimated crude Xome size"] = str(lc.calc_xome_size(throughput))

    with open(json_path, "w") as f:
        logger.info("Quality measurements were written into a JSON file: %s" % json_path)
        json.dump(tobe_json, f, indent=4)

    logger.info("Generated a json summary.")

    root_dict = {}

    root_dict['stats']  = OrderedDict()
    if suffix == "":
        root_dict['stats']['Sample name'] = "-"
    else:
        root_dict['stats']['Sample name'] = suffix.replace('_', '')
    root_dict['stats']['Yield'] = int(throughput)
    root_dict['stats']['Number of reads'] = len(lengths)

    if args.sequel or file_format_code == 3: # fasta has no qual
        root_dict['stats']['Q7 bases'] = "-"
    else:
        root_dict['stats']['Q7 bases'] = "%.3f%%" % float(100*q7/throughput)
    root_dict['stats']['Longest read'] = int(longest)

    if lc.get_unmapped_med_frac():
        root_dict['stats']['Estimated non-sense read fraction'] = "%.3f" % float(lc.get_unmapped_med_frac())

    if lc.get_control_frac():
        root_dict['stats']['Estimated spiked-in control read fraction'] = "%.3f" % float(lc.get_control_frac())

    if (args.adp5 and max_iden_adp5 >= 0.75) or (args.adp3 and max_iden_adp3 >= 0.75):
        root_dict['ad'] = OrderedDict()
    if args.adp5 and max_iden_adp5 >= 0.75:
        root_dict['ad']["Number of trimmed reads in 5\' "] = num_trim5
        root_dict['ad']["Max seq identity for the adpter in 5\'"] = "%.3f" % max_iden_adp5
        root_dict['ad']["Average trimmed length in 5\'"] = "%.3f" % np.mean(adp_pos5)
    if args.adp3 and max_iden_adp3 >= 0.75:
        root_dict['ad']["Number of trimmed reads in 3\'"] = num_trim3
        root_dict['ad']["Max seq identity for the adpter in 3\'"] = "%.3f" % max_iden_adp3
        root_dict['ad']["Average trimmed length in 3\'"] = "%.3f" % np.mean(adp_pos3)

    if args.pb:
        root_dict['pb'] = True

    if args.sequel :
        root_dict['sequel'] = True

    root_dict['rl'] = {'name':enc_b64_str(fig_path),\
                      'stats':OrderedDict([\
                               ('Mean read length', "%.3f" % mean_len),\
                               ('N50', "%.3f" % n50)])}
    root_dict['rq'] = {'name':enc_b64_str(fig_path_rq)}

    if lc.is_no_coverage():
        root_dict['rc'] = {'cov_plot_name':enc_b64_str(fig_path_cv),
                           'cov_over_len_plot_name':enc_b64_str(fig_path_cl),\
                           'cov_ovlp_qv_plot_name':enc_b64_str(fig_path_qv),\
                           'stats':OrderedDict([\
                                                ('Number of sampled reads', s_n_seqs),\
                                                ('Mean per read coverage', "N/A"),\
                                                ('S.D. per read coverage', "N/A"), \
                                                ('Crude estimated Xome size', lc.calc_xome_size(throughput))])}
    elif args.transcript:
        root_dict['rc'] = {'cov_plot_name':enc_b64_str(fig_path_cv),
                           'cov_over_len_plot_name':enc_b64_str(fig_path_cl),\
                           'cov_ovlp_qv_plot_name':enc_b64_str(fig_path_qv),\
                           'stats':OrderedDict([\
                                                ('Number of sampled reads', s_n_seqs),\
                                                ('Mode of per read coverage', "%.3f" % lc.get_logn_mode()),\
                                                ('mu of per read coverage', "%.3f" % lc.get_logn_mu()), \
                                                ('sigma of per read coverage', "%.3f" % lc.get_logn_sigma()), \
                                                ('Crude estimated Xome size', lc.calc_xome_size(throughput))])}

    elif lc.is_low_coverage():
        root_dict['rc'] = {'cov_plot_name':enc_b64_str(fig_path_cv),
                           'cov_over_len_plot_name':enc_b64_str(fig_path_cl),\
                           'cov_ovlp_qv_plot_name':enc_b64_str(fig_path_qv),\
                           'stats':OrderedDict([\
                                                ('Number of sampled reads', s_n_seqs),\
                                                ('Mode of per read coverage', "%.3f" % lc.get_logn_mode()),\
                                                ('mu of per read coverage', "%.3f" % lc.get_logn_mu()), \
                                                ('sigma of per read coverage', "%.3f" % lc.get_logn_sigma()), \
                                                ('Crude estimated Xome size', lc.calc_xome_size(throughput))])}

    else:
        root_dict['rc'] = {'cov_plot_name':enc_b64_str(fig_path_cv),
                           'cov_over_len_plot_name':enc_b64_str(fig_path_cl),\
                           'cov_ovlp_qv_plot_name':enc_b64_str(fig_path_qv),\
                           'stats':OrderedDict([\
                                                ('Number of sampled reads', s_n_seqs),\
                                                ('Mean per read coverage', "%.3f" % lc.get_mean()),\
                                                ('S.D. per read coverage', "%.3f" % lc.get_sd()), \
                                                ('Crude estimated Xome size', lc.calc_xome_size(throughput))])}

    root_dict['gc'] = {'name':enc_b64_str(fig_path_gc),\
                      'stats':OrderedDict([\
                               ('Mean per read GC content', "%.3f %%" % (100.0 * gc_read_mean)),\
                               ('s.d. per read GC content', "%.3f %%" % (100.0 * gc_read_sd))
                                       ])}
    root_dict['fr'] = {'name':enc_b64_str(fig_path_ta)}
    root_dict['sc'] = {'name':enc_b64_str(fig_path_ma)}

    # alerts
    root_dict['warns'] = OrderedDict()
    root_dict['errors'] = OrderedDict()

    if not args.sequel and file_format_code == 2: # pb-bam and fasta have no qual
        if q7/throughput <= 0.65 and q7/throughput > 0.5:
            root_dict['warns']['Low Q7'] = 'This value should be higher than 65%.'
        elif q7/throughput <= 0.5:
            root_dict['errors']['Too low Q7'] = 'This value should be higher than 50%. Ideally, higher than 65%.'

    if lc.is_no_coverage():
        root_dict['errors']['Coverage estimation failure'] = 'Coverage estimation cannot be made. No or very little coverage data exists.'
    elif very_low_coverage_mode:
        if lc.is_low_coverage():
            root_dict['warns']['Low coverage'] = 'Coverage of data looks to be very low/skewed.'
        else:
            root_dict['warns']['Low coverage'] = 'Coverage of data looks to be very low.'
        e_zero = lc.get_expected_zero_rate()
        logger.info("Low coverage mode: expected zero rate for the given coverage %.2f is %.2f." % e_zero )
        #adj_e = float(1.0 - len(mapped_ids)/float(s_n_seqs-lc.get_control_num()))
        adj_e = lc.get_unmapped_med_frac() - e_zero[1]
        if adj_e >= nonsense_read_warn_threshold and adj_e < nonsense_read_error_threshold:
            root_dict['warns']['High non-sense read fraction'] = 'This value should be lower than %.2f%%.' % float((nonsense_read_warn_threshold+e_zero[1])*100)
        elif adj_e >= nonsense_read_error_threshold:
            root_dict['errors']['Too high non-sense read fraction'] = 'This value should not be higher than %.2f%%.' % float((nonsense_read_error_threshold+e_zero[1])*100)
    else:
        if lc.get_unmapped_med_frac() >= nonsense_read_warn_threshold and lc.get_unmapped_med_frac() < nonsense_read_error_threshold:
            root_dict['warns']['High non-sense read fraction'] = 'This value should be lower than %d%%.' % int(nonsense_read_warn_threshold*100)
        elif lc.get_unmapped_med_frac() >= nonsense_read_error_threshold:
            root_dict['errors']['Too high non-sense read fraction'] = 'This value should not be higher than %d%%.' % int(nonsense_read_error_threshold*100)

    if num_trim5 and not args.pb:
        if num_trim5/len(lengths) <= 0.3:
            root_dict['warns']['Low number of adapter hits in 5\''] = 'This value should be higher than 30% if adapter sequences were not removed.'

    if lc.get_errors():
        for e in lc.get_errors():
            root_dict['errors'][e[0]] = e[1]

    if lc.get_warnings():
        for w in lc.get_warnings():
            root_dict['warns'][w[0]] = w[1]

    template_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "web_summary")
    env = Environment(loader=FileSystemLoader(template_dir, encoding='utf8'))
    tpl = env.get_template('web_summary.tpl.html')
    html = tpl.render( root_dict )
    with open(html_path, "wb") as f:
        f.write(html.encode('utf-8'))
    #if not os.path.isdir(os.path.join(args.out, "css")):
    #    os.makedirs(os.path.join(args.out, "css"), exist_ok=True)
    #if not os.path.isdir(os.path.join(args.out, "vendor")):
    #    os.makedirs(os.path.join(args.out, "vendor"), exist_ok=True)
    #if not os.path.isdir(os.path.join(args.out, "figs")):
    #    os.makedirs(os.path.join(args.out, "figs"), exist_ok=True)
    #copytree(os.path.join(template_dir, 'css'), os.path.join(args.out, "css"))
    #copytree(os.path.join(template_dir, 'vendor'), os.path.join(args.out, "vendor"))
    #copytree(os.path.join(template_dir, 'figs'), os.path.join(args.out, "figs"))
    logger.info("Generated a summary html.")

    logger.info("Finished all processes.")

    if args.db and file_format_code != 0:
        #tempdb_path = os.path.join(args.out, "analysis", "minimap2", "t_db_minimap2")
        if os.path.exists(tempdb_path):
            try:
                os.remove(tempdb_path)
                logger.info("tmpdb %s was removed." % tempdb_path)
            except (OSError, e):
                logger.error("%s - %s." % (e.filename, e.strerror))
        else:
            logger.warning("tmpdb file %s does not exist. skip removal of this file." % tempdb_path)
        if args.short:
            #tempdb_short_path = os.path.join(args.out, "analysis", "minimap2", "t_db_minimap2_short")
            if os.path.exists(tempdb_short_path):
                try:
                    os.remove(tempdb_short_path)
                    logger.info("tmpdb %s was removed." % tempdb_short_path)
                except (OSError, e):
                    logger.error("%s - %s." % (e.filename, e.strerror))
            else:
                logger.warning("tmpdb file %s does not exist. skip removal of this file." % tempdb_short_path)
        logger.info("Cleaned up temp db(s).")

# stand alone
if __name__ == "__main__":
    # parsing
    parser = argparse.ArgumentParser(
        prog='LongQC',
        description='LongQC is a software to asses the quality of long read data from the third generation sequencers.',
        add_help=True,
    )
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)
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
    presets = ["pb-rs2", "pb-sequel", "pb-hifi", "ont-ligation", "ont-rapid", "ont-1dsq"]
    help_preset = 'a platform/kit to be evaluated. adapter and some ovlp parameters are automatically applied. ('+", ".join(presets)+')'
    parser_sample = subparsers.add_parser('sampleqc', help='see `sampleqc -h`')
    parser_sample.add_argument('input', help='Input [fasta, fastq or pbbam]', type=str)
    parser_sample.add_argument('-o', '--output', \
                               help='path for output directory', type=str,\
                               dest = 'out', required=True, default = None)
    parser_sample.add_argument('-x', '--preset', choices=presets, help=help_preset, metavar='preset', required=True)

    parser_sample.add_argument('-t', '--transcript', \
                               help='applies the preset for transcripts, RNA or cDNA sequences', \
                               dest = 'transcript', action = 'store_true', default = None)
    parser_sample.add_argument('-n', '--n_sample',\
                               help='the number of sequences for sampling. (>0 and <=10000) [Default is 5000].', type=int,\
                               dest = 'nsample', default = 5000)
    parser_sample.add_argument('-s', '--sample_name', \
                               help='sample name is added as a suffix for each output file.', type=str,\
                               dest = 'suf', default = None)

    parser_sample.add_argument('-c', '--trim_output', type=str,\
                               help='path for trimmed reads. If this is not given, trimmed reads won\'t be saved.',\
                               dest = 'trim', default = None)
    parser_sample.add_argument('--adapter_5', help='adapter sequence for 5\'.', dest = 'adp5', default = None)
    parser_sample.add_argument('--adapter_3', help='adapter sequence for 3\'.', dest = 'adp3', default = None)


    #parser_sample.add_argument('-a', '--accurate',\
    #                           help='this turns on the more sensitive setting. More accurate but slower.', action = 'store_true',\
    #                           dest = 'acc', default = None)
    parser_sample.add_argument('-f', '--fast',\
                               help='this turns off sensitive setting. Faster but less accurate.', action = 'store_true',\
                               dest = 'fast', default = None)
    parser_sample.add_argument('-p', '--ncpu', help='the number of cpus for LongQC analysis [Default is 4. >=4 is required.]', type=int, dest = 'ncpu', default = 4)
    parser_sample.add_argument('-d', '--db',\
                               help='make minimap2 db in parallel to other tasks.',\
                               dest = 'db', action = 'store_true', default = False)
    #parser_sample.add_argument('-l', '--list',\
    #                           help='give IDs of reads for analysis. This replaces subsampling step. (>0 and <=50000)',\
    #                           dest = 'list', default = None)
    parser_sample.add_argument('-m', '--mem', type=float,\
                               help='memory limit for chunking. Please specify in gigabytes (>0 and <=2). [Default is 0.5]', \
                               dest = 'mem', default = 0.5)
    parser_sample.add_argument('-i', '--index', type=str,\
                               help='Give index size for minimap2 (-I) in bp. Reduce when running on a small memory machine.'+
                               'Default is 4G.', dest = 'inds', default = '4G')

    parser_sample.add_argument('-b', '--short',\
                               help='this turns on the highly sensitive setting for very short and erroneous reads (<500bp).', action = 'store_true',\
                               dest = 'short', default = None)
    parser_sample.add_argument('--pb', help=argparse.SUPPRESS,\
                               dest = 'pb', action = 'store_true', default = None)
    parser_sample.add_argument('--sequel', help=argparse.SUPPRESS,\
                               dest = 'sequel', action = 'store_true', default = None)
    parser_sample.add_argument('--ont', help=argparse.SUPPRESS,\
                               dest = 'ont', action = 'store_true', default = None)
    parser_sample.set_defaults(handler=command_sample)

    # help
    parser_help = subparsers.add_parser('help', help='see `help -h`')
    parser_help.add_argument('command', help='')
    parser_help.set_defaults(handler=command_help)
    args = parser.parse_args()
        
    main(args)
