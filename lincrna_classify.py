# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrew.sczesnak@med.nyu.edu"
__date__ = "9/19/2011"
__version__ = 0.0

from optparse import OptionParser

parser = OptionParser(usage="%prog [options] <transcripts.gtf> <transcripts.tracking>",
                      version="%prog " + str(__version__))

parser.add_option("-o",
                  dest="output_dir",
                  metavar="[./lincrna_out]",
                  default="./lincrna_out",
                  help="write output files to this directory")
                  
parser.add_option("-p",
                  dest="num_threads",
                  type="int",
                  metavar="[1]",
                  default=1,
                  help="number of threads used during analysis")

# (1) Size selection
parser.add_option("--min-exons",
                  dest="min_exons",
                  type="int",
                  metavar="2",
                  default=2,
                  help="minimum number of exons")
                  
parser.add_option("--min-size",
                  dest="min_size",
                  type="int",
                  metavar="200",
                  default=200,
                  help="minimum transcript length")
                  
# (2) Minimal read coverage threshold
parser.add_option("--min-cov",
                  dest="min_cov",
                  type="int",
                  metavar="3",
                  default=3,
                  help="minimum transcript coverage")
                  
# (3) Filter of known non-lincRNA annotations
parser.add_option("-f",
                  dest="filter",
                  action="append",
                  metavar="filter.gtf",
                  help="remove transcripts with exons overlapping these GTFs")
                  
# (4) Positive coding potential threshold
parser.add_option("--min-csf",
                  dest="min_csf",
                  type="float",
                  metavar="10.0",
                  default=10.0,
                  help="minimum CSF score of transcript")
                  
# (5) Known protein domain filter
parser.add_option("--exclude-pfam",
                  dest="exclude_pfam",
                  action="store_true",
                  default=False,
                  help="remove transcripts containing Pfam domains")
                  
options, args = parser.parse_args()

