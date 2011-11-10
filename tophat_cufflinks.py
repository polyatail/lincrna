# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrew.sczesnak@med.nyu.edu"
__date__ = "10/26/2011"
__version__ = 0.0

from optparse import OptionParser, OptionGroup
import sys
import os
from wrappers import tophat
from wrappers import cufflinks
from wrappers import samtools
from wrappers import igvtools
import _common
import tempfile

def mkwig(input_bam, output_wig, filter_less_than = 1):
    st = samtools.SAMTools()
    pileup = st.pileup(input_bam)
    out_fp = open(output_wig, "w")

    this_chrom = None

    for line in pileup:
        line_split = line.strip().split()

        if int (line_split[3]) < filter_less_than: continue

        if line_split[0] != this_chrom:
            out_fp.write("variableStep chrom=%s\n" % (line_split[0],))

            this_chrom = line_split[0]

        out_fp.write("%s\t%s\n" % (line_split[1], line_split[3]))

    out_fp.flush()

def parse_options(arguments):
    global options, args

    parser = OptionParser(usage="%prog [options] <assembly> <tophat.params>",
                          version="%prog " + str(__version__))

    parser.add_option("-o",
                      dest="output_dir",
                      metavar="[./th_cl_out]",
                      default="./th_cl_out",
                      help="write output files to this directory")

    parser.add_option("-p",
                      dest="num_threads",
                      type="int",
                      metavar="[1]",
                      default=1,
                      help="number of threads used during analysis")

    parser.add_option("-r",
                      dest="reference",
                      metavar="reference.gtf",
                      default=None,
                      help="reference annotation (i.e. refGene) GTF file")

    tophat_opts = OptionGroup(parser, "TopHat Options")

    tophat_opts.add_option("--skip-tophat",
                           dest="skip_tophat",
                           action="store_true",
                           default=False,
                           help="skip execution of TopHat")

    tophat_opts.add_option("-c",
                           dest="num_procs",
                           type="int",
                           metavar="[1]",
                           default=1,
                           help="number of concurrent TopHat processes")


    parser.add_option_group(tophat_opts)

    cuff_opts = OptionGroup(parser, "Cufflinks Options")

    cuff_opts.add_option("--skip-cufflinks",
                         dest="skip_cufflinks",
                         action="store_true",
                         default=False,
                         help="skip execution of Cufflinks")


    cuff_opts.add_option("-M",
                         dest="mask",
                         metavar="mask.gtf",
                         default=None,
                         help="ignore all alignment within transcripts in this file")

    parser.add_option_group(cuff_opts)

    tdf_opts = OptionGroup(parser, "TDF Options")

    tdf_opts.add_option("--skip-tdf",
                        dest="skip_tdf",
                        action="store_true",
                        default=False,
                        help="do not make TDF files for each sample")

    tdf_opts.add_option("--filter-less-than",
                        dest="filter_less_than",
                        action="store",
                        type="int",
                        metavar="[1]",
                        default=1,
                        help="Filter positions with less than this coverage")

    parser.add_option_group(tdf_opts)

    options, args = parser.parse_args(arguments)

    if len(args) <> 2:
        print "Error: Incorrect number of arguments"
        parser.print_help()
        sys.exit(0)

    if not (options.skip_tophat and options.skip_cufflinks) and not options.reference:
        print "Error: Reference GTF required"
        parser.print_help()
        sys.exit(0)

    labels = []
    runs = []

    for line in open(args[1], "r"):
        if line.startswith("#"): continue
        line_split = line.strip().split("\t")

        if len(line_split) <> 6:
            print "Error: Invalid line in params file"
            print "\t", line
            parser.print_help()
            sys.exit(0)

        labels.append(line_split[0])
        runs.append(dict(zip(["label", "inner_dist", "inner_dist_sd",
                              "seed_len", "left_reads", "right_reads"],
                             line_split)))

    options.labels = labels
    options.runs = runs
    options.bowtie_index = _common.bowtie_index[args[0]]

def main():
    parse_options(sys.argv[1:])

    if not os.path.exists(options.output_dir):
        os.mkdir(options.output_dir)

    # tophat runs
    if not options.skip_tophat:
        jobs = []

        for label, run_params in zip(options.labels, options.runs):
            run_out_dir = os.path.join(options.output_dir, label)
            run_input_fastq = [run_params["left_reads"]]
            extra_params = ["--seed-length", run_params["seed_len"]]

            if run_params["right_reads"] != "None":
                run_input_fastq.append(run_params["right_reads"])
                extra_params.extend(["--mate-inner-dist", run_params["inner_dist"],
                                     "--mate-std-dev", run_params["inner_dist_sd"]])

            jobs.append((run_input_fastq,
                         ["--GTF", options.reference,
                          "--no-novel-juncs",
                          "--no-novel-indels",
                          "--no-coverage-search"] + \
                         extra_params,
                         run_out_dir))

        th = tophat.TopHat(options.bowtie_index,
                           options.output_dir,
                           options.num_threads)
        th.run_multiple(jobs, options.num_procs)

    # cufflinks runs
    if not options.skip_cufflinks:
        for label, run_params in zip(options.labels, options.runs):
            input_bam = os.path.join(options.output_dir, label, "accepted_hits.bam")

            if not os.path.exists(input_bam):
                print "Error: Could not find %s" % (input_bam,)
                continue

            run_out_dir = os.path.join(options.output_dir, "FPKM_" + label)

            cl = cufflinks.Cufflinks(run_out_dir, options.num_threads)
            cl.run(input_bam, ["--min-frags-per-transfrag", "0",
                               "-G", options.reference] + \
                              ([] if not options.mask else ["-M", options.mask]))

    # make TDF files
    if not options.skip_tdf:
        # collect sample data
        sample_info = []

        for label, run_params in zip(options.labels, options.runs):
            input_bam = os.path.join(options.output_dir, label, "accepted_hits.bam")

            if not os.path.exists(input_bam):
                print "Error: Could not find %s" % (input_bam,)
                continue

            wig_file = tempfile.NamedTemporaryFile("w", suffix=".wig", delete=False)
            tdf_out = os.path.join(options.output_dir, label, label + ".tdf")

            sample_info.append((input_bam, wig_file, tdf_out))

        # make wigs
        st = samtools.SAMTools()
        pid_list = []

        for input_bam, wig_file, tdf_out in sample_info:
            pid = _common.fork_and_run(st.wig,
                                       input_bam,
                                       wig_file,
                                       options.filter_less_than)

            pid_list.append(pid)

            _common.wait_for_slot(pid_list, options.num_threads)

        _common.wait_for_slot(pid_list, 0, True)

        # make tdfs
        igv = igvtools.IGVTools()
        pid_list = []
 
        for input_bam, wig_file, tdf_out in sample_info:
            pid = _common.fork_and_run(igv.tile,
                                       wig_file.name,
                                       tdf_out,
                                       args[0])

            pid_list.append(pid)

            _common.wait_for_slot(pid_list, options.num_threads)

        _common.wait_for_slot(pid_list, 0, True)

        # remove tempfiles
        for _, wig_file, _ in sample_info:
            os.unlink(wig_file.name)

if __name__ == "__main__":
    main()
