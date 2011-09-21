# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrew.sczesnak@med.nyu.edu"
__date__ = "9/19/2011"
__version__ = 0.0

from optparse import OptionParser
import copy
import shlex
import pprint

def _gtf_parser(gtf_file):
    empty_tx = {"raw": {},
                "exonStarts": {},
                "exonEnds": {},
                "chrom": None,
                "tx_len": 0,
                "num_exons": None,
                "csf": None,
                "cov": None,
                "pfams": None}
                
    transcripts = {}
    
    with open(gtf_file, "r") as gtf_fp:
        for line_num, line in enumerate(gtf_fp):
            line_split = line.strip().split("\t")
            
            metadata = shlex.split(line_split[8])
            metadata = dict([(metadata[i], metadata[i+1].rstrip(";")) for i in \
                       range(0, len(metadata)-1, 2)])
                       
            tx_id = metadata["transcript_id"]
            exon_num = metadata["exon_number"]
            
            try:
                transcripts[tx_id]
            except KeyError:
                transcripts[tx_id] = copy.deepcopy(empty_tx)

            transcripts[tx_id]["raw"][exon_num] = line            
            transcripts[tx_id]["exonStarts"][exon_num] = int(line_split[3])
            transcripts[tx_id]["exonEnds"][exon_num] = int(line_split[4])

            if transcripts[tx_id]["chrom"] == None:
                transcripts[tx_id]["chrom"] = line_split[0]
            elif transcripts[tx_id]["chrom"] != line_split[0]:
                raise ValueError("Chromosome mismatch in %s, line %s" % (tx_id, line_num))
            
    for tx_id, tx_data in transcripts:
        if len(tx_data["exonStarts"]) <> len(tx_data["exonEnds"]):
            raise ValueError("Exon count mismatch in %s" % (tx_id,))
            
        transcripts[tx_id]["num_exons"] = len(tx_data["exonStarts"])
        
        for start, end in zip(tx_data["exonStarts"], tx_data["exonEnds"]):
            if start > end:
                raise ValueError("Exon start > end in %s" % (tx_id,))
                
            transcripts[tx_id]["tx_len"] += end - start
            
    return transcripts
            
def main():
    _gtf_parser(args[0])
            
if __name__ == "__main__":
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

    main()