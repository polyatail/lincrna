# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrew.sczesnak@med.nyu.edu"
__date__ = "9/20/2011"
__version__ = 0.0

from optparse import OptionParser
import sys
import os
import _common
from wrappers import cufflinks, cuffcompare, samtools, scripture

def parse_options(arguments):
    global options, args

    parser = OptionParser(usage="%prog [options] <assembly> <sample1.sam> [...sampleN.sam]",
                          version="%prog " + str(__version__))
                          
    parser.add_option("-o",
                      dest="output_dir",
                      metavar="[./txome_out]",
                      default="./txome_out",
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
                      dest="reference",
                      metavar="reference.gtf",
                      default=None,
                      help="reference annotation (i.e. refGene) GTF file")
                      
    parser.add_option("-M",
                      dest="mask",
                      metavar="mask.gtf",
                      default=None,
                      help="ignore all alignment within transcripts in this file")

    parser.add_option("--prerun",
                      dest="prerun",
                      action="store_true",
                      default=False,
                      help="perform only the 'discovery' Cufflinks/Scripture runs")

    parser.add_option("--pool-transcripts",
                      dest="pool_transcripts",
                      action="store_true",
                      default=False,
                      help="pool transcripts with Cuffcompare")

    parser.add_option("--realrun",
                      dest="realrun",
                      action="store_true",
                      default=False,
                      help="perform only the quantitation Cufflinks runs")
    
    options, args = parser.parse_args()
    
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

    options.genome_dir = _common.genome_fasta[args[0]]    
    options.assembly, options.chroms = _common.find_genome_files(options.genome_dir)

    if options.assembly != args[0]:
        print "Error: Specified assembly does not match assembly found"
        parser.print_help()
        sys.exit(0)

def main():
    parse_options(sys.argv[1:])

    normal_run = not (options.realrun or options.prerun or options.pool_transcripts)

    if not os.path.exists(options.output_dir):
        os.mkdir(options.output_dir)

    if normal_run or options.prerun:
        # run cufflinks in 'discovery' mode
        for label, input_bam in zip(options.labels, args[1:]):
            run_out_dir = os.path.join(options.output_dir, "cufflinks_" + label)
                                      
            cl = cufflinks.Cufflinks(run_out_dir, options.num_threads)
            cl.run(input_bam, ["--min-frags-per-transfrag", "0"] + \
                              ([] if not options.mask else ["-M", options.mask]))
            
        # make sure there's a BAM index on every input file
        st = samtools.SAMTools()

        for input_bam in args[1:]:
            if st.index(input_bam) == False:
                raise ValueError("Failed to index %s" % (input_bam,))
            
        # run scripture in 'discovery' mode
        for label, input_bam in zip(options.labels, args[1:]):
            run_out_dir = os.path.join(options.output_dir, "scripture_" + label)
            run_out_bed = os.path.join(run_out_dir, "transcripts.bed")
            run_out_gtf = os.path.join(run_out_dir, "transcripts.gtf")
            
            sc = scripture.Scripture(options.genome_dir,
                                     run_out_dir,
                                     options.num_threads)
            sc.all_chrom(input_bam)
            scripture.bed_to_gtf(open(run_out_bed, "r"),
                                 open(run_out_gtf, "w"))

    # run cuffcompare to generate the unique intersection of all output GTFs
    if normal_run or options.pool_transcripts:
        input_gtfs = []
        
        for label in options.labels:
            input_gtfs.append(os.path.join(options.output_dir,
                                           "cufflinks_" + label,
                                           "transcripts.gtf"))
            input_gtfs.append(os.path.join(options.output_dir,
                                           "scripture_" + label,
                                           "transcripts.gtf"))

        run_out_dir = os.path.join(options.output_dir, "cuffcompare")
        
        cc = cuffcompare.Cuffcompare(options.genome_dir,
                                     run_out_dir)
        cc.run(input_gtfs,
               ([] if not options.reference else ["-r", options.reference]))
    
    # re-run each sample separately in abundance calculation mode
    if normal_run or options.realrun:
        for label, input_bam in zip(options.labels, args[1:]):
            run_out_dir = os.path.join(options.output_dir, "FPKM_" + label)
            
            cl = cufflinks.Cufflinks(run_out_dir, options.num_threads)
            cl.run(input_bam, ["--min-frags-per-transfrag", "0",
                               "-G", os.path.join(options.output_dir,
                                                  "cuffcompare",
                                                  "cuffcompare.combined.gtf")] + \
                              ([] if not options.mask else ["-M", options.mask]))
                        
if __name__ == "__main__":
    main()
