# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrew.sczesnak@med.nyu.edu"
__date__ = "9/19/2011"
__version__ = 0.0

from optparse import OptionParser
from Bio import SeqIO
import sys
import os

def _illumina_version(fname):
    first_record = SeqIO.parse(fname, "fastq").next()
    desc_split = first_record.description.split(" ")
    
    if len(desc_split) == 2:
        # Illumina v1.8
        illumina_ver = _illumina18
        
        meta_split = desc_split[1].split(":")
        
        if len(meta_split) <> 4:
            raise ValueError("Illumina v1.8+ metadata format invalid")
    elif len(desc_split) == 1:
        # Illumina v1.4
        illumina_ver = _illumina14
        
        meta_split = desc_split[0].split(":")
        
        if len(meta_split) <> 8:
            raise ValueError("Illumina v1.4 metadata format invalid")
    else:
        raise ValueError("Could not detect Illumina pipeline version")
        
    return illumina_ver

def _illumina14(description):
    meta_split = description.split(":")

    try:
        # IMPORTANT: 1 = KEEP THE READ, 0 = DISCARD THE READ
        # INTERNALLY, Y = DISCARD THE READ, N = KEEP THE READ
        if meta_split[7] == "1":
            meta_split[7] = "N"
        elif meta_split[7] == "0":
            meta_split[7] = "Y"
        else:
            raise ValueError("Filtered field must be 1/0")
    except IndexError:
        raise ValueError("Illumina v1.4 metadata format invalid")
    
    return meta_split[6], meta_split[7], ":".join(meta_split[:6])
    
def _illumina18(description):
    try:
        main_split = description.split(" ")
        meta_split = main_split[1].split(":")
    except IndexError:
        raise ValueError("Illumina v1.8+ metadata format invalid")
        
    if len(meta_split) <> 4:
        raise ValueError("Illumina v1.8+ metadata format invalid")
    
    return meta_split[0], meta_split[1], main_split[0]

def _single_parser(fp_in, fp_out, callback):
    for seq_rec in SeqIO.parse(fp_in, "fastq"):
        mate_pair, filtered, _ = callback(seq_rec.description)
        
        if mate_pair != "1":
            raise ValueError("Found mate_pair = 2 in single-end library")
        
        if filtered == "N":
            SeqIO.write(seq_rec, fp_out, "fastq")
        elif filtered == "Y":
            pass
        else:
            raise ValueError("Filtered field must be Y/N")
            
def _paired_parser(fp_in, fp_out_left, fp_out_right, callback):      
    left_count = 0
    right_count = 0
    
    for seq_rec in SeqIO.parse(fp_in, "fastq"):
        mate_pair, filtered, _ = callback(seq_rec.description)

        if filtered == "N":
            if mate_pair == "1":
                left_count += 1
                SeqIO.write(seq_rec, fp_out_left, "fastq")
            elif mate_pair == "2":
                right_count += 1
                SeqIO.write(seq_rec, fp_out_right, "fastq")
            else:
                raise ValueError("Paired end field must be 1/2")
        elif filtered == "Y":
            pass
        else:
            raise ValueError("Filtered field must be Y/N")
            
    fp_out_left.seek(0)
    fp_out_right.seek(0)
    
    left_parser = SeqIO.parse(fp_out_left, "fastq")
    right_parser = SeqIO.parse(fp_out_right, "fastq")

    if left_count <> right_count:
        raise ValueError("Left read count (%s) != right read count (%s)" % \
                         (left_count, right_count))

    for _ in range(left_count):
        left_rec = left_parser.next()
        right_rec = right_parser.next()

        if callback(left_rec.description)[2] != callback(right_rec.description)[2]:
            raise ValueError("Output reads are not in order")
        
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
        raise ValueError("Incorrect number of arguments")

def main():
    parse_options(sys.argv[1:])

    illumina_ver = _illumina_version(args[0])
    
    if not os.path.exists(options.output_dir):
        os.mkdir(options.output_dir)

    stripped_fname = ".".join(args[0].split(".")[:-1])

    if options.paired_end:
        left_out = os.path.join(options.output_dir,
                                stripped_fname + "-left.fastq")
        right_out = os.path.join(options.output_dir,
                                 stripped_fname + "-right.fastq")
                            
        _paired_parser(open(args[0], "r"),
                       open(left_out, "w+"),
                       open(right_out, "w+"),
                       illumina_ver)
    else:
        filtered_out = os.path.join(options.output_dir,
                                    stripped_fname + "-filtered.fastq")

        _single_parser(open(args[0], "r"),
                       open(filtered_out, "w+"),
                       illumina_ver)
        
if __name__ == "__main__":
    main()