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
    first_record = SeqIO.parse(args[0], "fastq").next()
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
    else:
        raise ValueError("Could not detect Illumina pipeline version")
        
    return illumina_ver

def _illumina14(description):
    meta_split = description.split(":")
    
    if meta_split[7] == "1":
        meta_split[7] = "Y"
    elif meta_split[7] == "0":
        meta_split[7] = "N"
    else:
        raise ValueError("Filtered field must be 1/0")
    
    return meta_split[6], meta_split[7], ":".join(meta_split[:6])
    
def _illumina18(description):
    main_split = description.split(" ")
    meta_split = main_split[1].split(":")
    
    return meta_split[0], meta_split[1], main_split[0]

def _single_parser(fname, callback):
    stripped_fname = ".".join(fname.split(".")[:-1])
    filtered_out = open(os.path.join(options.output_dir,
                                     stripped_fname + "-filtered.fastq"), "w")

    for seq_rec in SeqIO.parse(args[0], "fastq"):
        mate_pair, filtered, _ = callback(seq_rec.description)
        
        assert mate_pair == "1"
        
        if filtered == "Y":
            SeqIO.write(seq_rec, filtered_out, "fastq")
        elif filtered == "N":
            pass
        else:
            raise ValueError("Filtered field must be Y/N")
            
def _paired_parser(fname, callback):
    stripped_fname = ".".join(fname.split(".")[:-1])
    left_out = open(os.path.join(options.output_dir,
                                 stripped_fname + "-left.fastq"), "w+")
    right_out = open(os.path.join(options.output_dir,
                                  stripped_fname + "-right.fastq"), "w+")
                                  
    left_count = 0
    right_count = 0
    
    for seq_rec in SeqIO.parse(args[0], "fastq"):
        mate_pair, filtered, _ = callback(seq_rec.description)

        if filtered == "Y":
            if mate_pair == "1":
                left_count += 1
                SeqIO.write(seq_rec, left_out, "fastq")
            elif mate_pair == "2":
                right_count += 1
                SeqIO.write(seq_rec, right_out, "fastq")
            else:
                raise ValueError("Paired end field must be 1/2")
        elif filtered == "N":
            pass
        else:
            raise ValueError("Filtered field must be Y/N")
            
    left_out.seek(0)
    right_out.seek(0)
    
    left_parser = SeqIO.parse(left_out, "fastq")
    right_parser = SeqIO.parse(right_out, "fastq")

    if left_count <> right_count:
        raise ValueError("Left read count (%s) != right read count (%s)" % \
                         (left_count, right_count))

    for _ in range(left_count):
        left_rec = left_parser.next()
        right_rec = right_parser.next()

        assert callback(left_rec.description)[2] == callback(right_rec.description)[2]

def main():
    illumina_ver = _illumina_version(args[0])
    
    if not os.path.exists(options.output_dir):
        os.mkdir(options.output_dir)

    if options.paired_end:
        _paired_parser(args[0], illumina_ver)
    else:
        _single_parser(args[0], illumina_ver)
        
if __name__ == "__main__":
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
              
    options, args = parser.parse_args()
    
    if len(args) <> 1:
        raise ValueError("Incorrect number of arguments")

    main()