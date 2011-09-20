# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrew.sczesnak@med.nyu.edu"
__date__ = "9/19/2011"
__version__ = 0.0

from optparse import OptionParser
import os
import subprocess
import time

def run_tophat(prefix, tophat_options, sample_reads, sample_name):
    sample_reads = sample_reads.split(",")
    outdir = os.path.join(options.output_dir, prefix + sample_name)

    if len(sample_reads) == 2:
        paired_end_args = ["-r", str(options.inner_dist),
                           "--mate-std-dev", str(options.inner_dist_sd)]
    elif len(cond_reads) > 2:
        raise ValueError("Maximum of two reads per sample!")
    else:
        paired_end_args = []

    if os.path.exists(outdir):
        raise ValueError("Output directory %s already exists!" % (outdir,))

    os.mkdir(outdir)
    
    with open(os.path.join(outdir, "tophat.log"), "w") as th_log:
        tophat_proc = subprocess.Popen(["tophat",
                                        "-p", str(options.num_threads),
                                        "-o", outdir,
                                        "-z", "none"] + \
                                        tophat_options + \
                                        paired_end_args + \
                                        [args[0]] + \
                                        sample_reads,
                                        stderr=th_log)
                                        
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
            if line.startswith("track"): continue

            line_split = line.split("\t")

            juncs.append("\t".join(line_split[:3] + [line_split[5]]))
            
    return juncs

def main():
    if not os.path.exists(options.output_dir):
        os.mkdir(options.output_dir)
    
    # run tophat on each sample individually, using default params with
    # --min-anchor=5 and --min-isoform-fraction=0
    for cond_name, cond_reads in zip(options.labels, args[1:]):
        run_tophat("prerun_", ["-F", "0", "-a", "5"], cond_reads, cond_name)

    # generate pooled junctions across all samples
    pooled_juncs_file = os.path.join(options.output_dir, "pooled.juncs")
    pooled_juncs = []
    
    for cond_name in options.labels:
        pooled_juncs.extend(_bed_to_junc(os.path.join(options.output_dir,
                                                      "prerun_" + sample_name,
                                                      "junctions.bed")))

    with open(pooled_juncs_file, "w") as fp:
        fp.write("\n".join(pooled_juncs))
        
    # re-run tophat on each sample individually, using default params with
    # --raw-juncs and --no-novel-juncs
    for cond_name, cond_reads in zip(options.labels, args[1:]):
        run_tophat("", ["-j", pooled_juncs_file,
                        "--no-novel-juncs"], cond_reads, cond_name)
                        
if __name__ == "__main__":
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
    
    options, args = parser.parse_args()
    
    if len(args) < 2:
        raise ValueError("Not enough arguments")
        
    if options.labels != None:
        options.labels = options.labels.split(",")
        
        if len(options.labels) <> len(args) - 1:
            raise ValueError("When using -L, must specify a label for every condition")
    else:
        options.labels = ["sample%s" % (x,) for x in range(len(args))]
        
    main()