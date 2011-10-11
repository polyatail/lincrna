# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrew.sczesnak@med.nyu.edu"
__date__ = "9/19/2011"
__version__ = 0.0

from optparse import OptionParser
import sys
import os
from lib import fastq
        
def parse_options(arguments):
    global options, args

    parser = OptionParser(usage="%prog [options] <reads.fastq>",
                          version="%prog " + str(__version__))
                          
    parser.add_option("-o",
                      dest="output_dir",
                      metavar="[./filter_out]",
                      default="./filter_out",
                      help="write output files to this directory")
                      
    parser.add_option("--paired-end",
                      dest="paired_end",
                      action="store_true",
                      default=False,
                      help="this is a paired-end library")
              
    options, args = parser.parse_args(arguments)
    
    if len(args) <> 1:
        print "Error: Incorrect number of arguments"
        parser.print_help()
        sys.exit(0)

def main():
    parse_options(sys.argv[1:])
    
    fastq_ver = fastq.fastq_version(args[0])
    fastq_readlen = fastq.fastq_readlen(args[0])

    callback_func = fastq.fastq_ver_to_callback(fastq_ver)
    phred_offset = fastq.fastq_ver_to_phred(fastq_ver)
    
    if not os.path.exists(options.output_dir):
        os.mkdir(options.output_dir)

    stripped_fname = ".".join(os.path.basename(args[0]).split(".")[:-1])

    if options.paired_end:
        left_out = os.path.join(options.output_dir,
                                stripped_fname + "-left.fastq")
        right_out = os.path.join(options.output_dir,
                                 stripped_fname + "-right.fastq")
        orphan_out = os.path.join(options.output_dir,
                                  stripped_fname + "-orphans.fastq")
                            
        fastq.paired_parser(open(args[0], "r"),
                            open(left_out, "w+"),
                            open(right_out, "w+"),
                            callback_func,
                            phred_offset,
                            25,
                            int(fastq_readlen) / 2)
    else:
        filtered_out = os.path.join(options.output_dir,
                                    stripped_fname + "-filtered.fastq")

        fastq.single_parser(open(args[0], "r"),
                            open(filtered_out, "w+"),
                            callback_func,
                            phred_offset,
                            25,
                            int(fastq_readlen) / 2)
        
if __name__ == "__main__":
    main()
