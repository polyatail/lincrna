# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrew.sczesnak@med.nyu.edu"
__date__ = "9/19/2011"
__version__ = 0.0

from optparse import OptionParser
import _common
import filter_reads
import sys
import os
import subprocess
import time

def _validate_reads(sample_reads):
    # figure out phred type and read length for each file
    fastq_ver = []
    phred_ver = []
    readlen = []

    for fname in sample_reads:
        fq_ver = filter_reads.fastq_version(fname)
        ph_ver = filter_reads.phred_version(fq_ver)
        rl = filter_reads.fastq_readlen(fname)

        fastq_ver.append(fq_ver)
        phred_ver.append(ph_ver)
        readlen.append(rl)

    # make sure they're consistent
    if len(set(fastq_ver)) <> 1:
        raise ValueError("Passed reads files have multiple FASTQ types: %s" % \
                         (fastq_ver,))
    else:
        fastq_ver = fastq_ver[0]

    if len(set(phred_ver)) <> 1:
        raise ValueError("Passed reads files have multiple phred types: %s" % \
                         (phred_ver,))
    else:
        phred_ver = phred_ver[0]

    if len(set(readlen)) <> 1:
        raise ValueError("Passed reads files have multiple read lengths: %s" % \
                         (readlen,))
    else:
        readlen = readlen[0]

    return fastq_ver, phred_ver, readlen

def run_tophat(prefix, tophat_options, sample_reads, sample_name):
    sample_reads = sample_reads.split(",")
    outdir = os.path.join(options.output_dir, prefix + sample_name)

    if len(sample_reads) == 2:
        paired_end_args = ["-r", str(options.inner_dist),
                           "--mate-std-dev", str(options.inner_dist_sd)]
    elif len(sample_reads) > 2:
        raise ValueError("Maximum of two reads per sample!")
    else:
        paired_end_args = []

    _, phred_ver, readlen = _validate_reads(sample_reads)

    # set segment length
    if readlen < 50:
        seg_length = int(readlen / 2)
    else:
        seg_length = 25

    if os.path.exists(outdir):
        raise ValueError("Output directory %s already exists!" % (outdir,))

    os.mkdir(outdir)
    
    with open(os.path.join(outdir, "tophat.log"), "w") as th_log:
        tophat_cmd = [_common.TOPHAT_PATH,
                      "-p", str(options.num_threads),
                      "-o", outdir,
                      "-z", "none",
                      "--segment-length", str(seg_length),
                      "--" + phred_ver] + \
                      tophat_options + \
                      paired_end_args + \
                      [args[0]] + \
                      sample_reads
                                        
        th_log.write(" ".join(tophat_cmd) + "\n\n")
        
        tophat_proc = subprocess.Popen(tophat_cmd, stderr=th_log)
                                        
    while tophat_proc.poll() == None:
        time.sleep(1)
        
    if os.path.exists(os.path.join(outdir, "junctions.bed")):
        return True
    else:
        return False
        
def _bed_to_junc(bedfile):
    juncs = []
    
    with open(bedfile, "r") as fp:
        for line in fp:
            if line.startswith("track") or not line.strip():
                continue

            line_split = line.split("\t")

            try:
                juncs.append("\t".join(line_split[:3] + [line_split[5]]))
            except IndexError:
                raise ValueError("Invalid BED file: %s" % (bedfile,))
            
    return juncs
    
def _pool_juncs(in_bedfiles, out_juncfile):
    pooled_juncs = []
    
    for bedfile in in_bedfiles:
        pooled_juncs.extend(_bed_to_junc(bedfile))

    with open(out_juncfile, "w") as fp:
        fp.write("\n".join(pooled_juncs))    
        fp.write("\n")

def parse_options(arguments):
    global options, args

    parser = OptionParser(usage="%prog [options] <bowtie_index> <reads1-L[,reads1-R]> [reads2-L[,reads2-R]]",
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
    
    options, args = parser.parse_args(arguments)
    
    if len(args) < 2:
        print "Error: Not enough arguments"
        parser.print_help()
        sys.exit(0)
        
    if options.labels != None:
        options.labels = options.labels.split(",")
        
        if len(options.labels) <> len(args) - 1:
            print "When using -L, must specify a label for every condition"
            parser.print_help()
            sys.exit(0)
    else:
        options.labels = ["sample%s" % (x,) for x in range(len(args))]

def main(arguments=sys.argv[1:]):
    parse_options(arguments)
    
    if not os.path.exists(options.output_dir):
        os.mkdir(options.output_dir)
    
    # run tophat on each sample individually, using default params with
    # --min-anchor=5 and --min-isoform-fraction=0
    for cond_name, cond_reads in zip(options.labels, args[1:]):
        run_tophat("prerun_",
                   ["-F", "0", "-a", "5"],
                   cond_reads,
                   cond_name)

    # generate pooled junctions across all samples
    _pool_juncs([os.path.join(options.output_dir,
                              "prerun_" + x,
                              "junctions.bed") for x in options.labels],
                os.path.join(options.output_dir, "pooled.juncs"))
        
    # re-run tophat on each sample individually, using default params with
    # --raw-juncs and --no-novel-juncs
    for cond_name, cond_reads in zip(options.labels, args[1:]):
        run_tophat("",
                   ["-j", pooled_juncs_file, "--no-novel-juncs"],
                   cond_reads,
                   cond_name)
                        
if __name__ == "__main__":
    main()
