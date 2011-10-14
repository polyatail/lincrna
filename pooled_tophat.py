# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrew.sczesnak@med.nyu.edu"
__date__ = "9/19/2011"
__version__ = 0.0

from optparse import OptionParser
import sys
import os
from wrappers import tophat
import _common

def pool_juncs(in_bedfiles, out_juncfile):
    pooled_juncs = []
    
    for bedfile in in_bedfiles:
        pooled_juncs.extend(tophat.bed_to_junc(bedfile))

    # make lines unique
    pooled_juncs = sorted(list(set(pooled_juncs)))

    with open(out_juncfile, "w") as fp:
        fp.write("\n".join(pooled_juncs))    
        fp.write("\n")

def parse_options(arguments):
    global options, args

    parser = OptionParser(usage="%prog [options] <assembly> <reads1-L[,reads1-R]> [reads2-L[,reads2-R]]",
                          version="%prog " + str(__version__))
                          
    parser.add_option("-o",
                      dest="output_dir",
                      metavar="[./tophat_out]",
                      default="./tophat_out",
                      help="write output files to this directory")
                      
    parser.add_option("-p",
                      dest="num_threads",
                      type="int",
                      metavar="[1]",
                      default=1,
                      help="number of threads used during analysis")

    parser.add_option("-c",
                      dest="num_procs",
                      type="int",
                      metavar="[1]",
                      default=1,
                      help="number of concurrent TopHat processes")
                      
    parser.add_option("-L",
                      dest="labels",
                      metavar="sample1,...sampleN",
                      default=None,
                      help="comma-separated list of condition labels")
                      
    parser.add_option("-r",
                      dest="inner_dist",
                      type="int",
                      metavar="[200]",
                      default=200,
                      help="the mean inner distance between mate-pairs")
                      
    parser.add_option("-s",
                      dest="inner_dist_sd",
                      type="int",
                      metavar="[20]",
                      default=20,
                      help="the standard deviation of the distance between pairs")

    parser.add_option("--prerun",
                      dest="prerun",
                      action="store_true",
                      default=False,
                      help="perform only the 'discovery' TopHat runs")

    parser.add_option("--pool-juncs",
                      dest="pool_juncs",
                      action="store_true",
                      default=False,
                      help="only pool junctions across samples")

    parser.add_option("--realrun",
                      dest="realrun",
                      action="store_true",
                      default=False,
                      help="perform only the alignment TopHat runs")
    
    options, args = parser.parse_args(arguments)
    
    if len(args) < 2:
        print "Error: Not enough arguments"
        parser.print_help()
        sys.exit(0)
        
    if options.labels != None:
        options.labels = options.labels.split(",")
        
        if len(options.labels) <> len(args) - 1:
            print "Error: When using -L, must specify a label for every condition"
            parser.print_help()
            sys.exit(0)
    else:
        options.labels = ["sample%s" % (x,) for x in range(len(args))]
        
    options.bowtie_index = _common.bowtie_index[args[0]]

def main(arguments=sys.argv[1:]):
    parse_options(arguments)
    
    if not os.path.exists(options.output_dir):
        os.mkdir(options.output_dir)

    normal_run = not (options.realrun or options.prerun or options.pool_juncs)
    
    # run tophat on each sample individually, using default params with
    # --min-anchor=5 and --min-isoform-fraction=0
    if normal_run or options.prerun:
        jobs = []

        for label, input_fastq in zip(options.labels, args[1:]):
            run_out_dir = os.path.join(options.output_dir, "prerun_" + label)
            run_input_fastq = input_fastq.split(",")
            
            jobs.append((run_input_fastq,
                         ["-F", "0", "-a", "5"] + \
                         ["-r", options.inner_dist,
                          "-s", options.inner_dist_sd],
                         run_out_dir))
            
        th = tophat.TopHat(options.bowtie_index,
                           options.output_dir,
                           options.num_threads)
        th.run_multiple(jobs, options.num_procs)

    # generate pooled junctions across all samples
    if normal_run or options.pool_juncs:
        pooled_juncs_file = os.path.join(options.output_dir, "pooled.juncs")
    
        pool_juncs([os.path.join(options.output_dir,
                                 "prerun_" + x,
                                 "junctions.bed") for x in options.labels],
                   pooled_juncs_file)
    else:
        pooled_juncs_file = os.path.join(options.output_dir, "pooled.juncs")
        
    # re-run tophat on each sample individually, using default params with
    # --raw-juncs and --no-novel-juncs
    if normal_run or options.realrun:
        jobs = []

        for label, input_fastq in zip(options.labels, args[1:]):
            run_out_dir = os.path.join(options.output_dir, label)
            run_input_fastq = input_fastq.split(",")
            
            jobs.append((run_input_fastq,
                         ["-j", pooled_juncs_file, "--no-novel-juncs"] + \
                         ["-r", options.inner_dist,
                          "-s", options.inner_dist_sd],
                         run_out_dir))
            
        th = tophat.TopHat(options.bowtie_index,
                           options.output_dir,
                           options.num_threads)
        th.run_multiple(jobs, options.num_procs)
                        
if __name__ == "__main__":
    main()
