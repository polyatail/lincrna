# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrew.sczesnak@med.nyu.edu"
__date__ = "9/19/2011"
__version__ = 0.0

from optparse import OptionParser
import sys
import copy
import shlex
import os
import _phylocsf

def _gtf_parser(gtf_file):
    empty_tx = {"raw": {},
                "exonStarts": {},
                "exonEnds": {},
                "exonClassCodes": {},
                "chrom": None,
                "num_exons": None,
                "tx_len": 0,    
                "tx_classcode": None,
                "annotation": None,
                "locus_id": None,
                "csf": None,
                "cov": {},
                "pfams": None,
                "note": None}
                
    transcripts = {}
    
    with open(gtf_file, "r") as gtf_fp:
        for line_num, line in enumerate(gtf_fp):
            line_split = line.strip().split("\t")
            
            if line_split[2] != "exon":
                continue
            
            metadata = shlex.split(line_split[8])
            metadata = dict([(metadata[i], metadata[i+1].rstrip(";")) for i in \
                       range(0, len(metadata)-1, 2)])
                       
            tx_id = metadata["transcript_id"]
            
            try:
                transcripts[tx_id]
            except KeyError:
                transcripts[tx_id] = copy.deepcopy(empty_tx)
                
            try:
                exon_num = metadata["exon_number"]
            except KeyError:
                exon_num = len(transcripts[tx_id]["exonStarts"]) + 1

            transcripts[tx_id]["raw"][exon_num] = line            
            transcripts[tx_id]["exonStarts"][exon_num] = int(line_split[3])
            transcripts[tx_id]["exonEnds"][exon_num] = int(line_split[4])
            
            try:
                transcripts[tx_id]["exonClassCodes"][exon_num] = metadata["class_code"]
            except KeyError:
                transcripts[tx_id]["exonClassCodes"][exon_num] = "."

            if transcripts[tx_id]["chrom"] == None:
                transcripts[tx_id]["chrom"] = line_split[0]
            elif transcripts[tx_id]["chrom"] != line_split[0]:
                raise ValueError("Chromosome mismatch in %s, line %s" % (tx_id, line_num))
            
    for tx_id, tx_data in transcripts.items():
        if len(tx_data["exonStarts"]) <> len(tx_data["exonEnds"]):
            raise ValueError("Exon count mismatch in %s" % (tx_id,))
            
        transcripts[tx_id]["num_exons"] = len(tx_data["exonStarts"])
        
        for start, end in zip(tx_data["exonStarts"].values(), tx_data["exonEnds"].values()):
            if start > end:
                raise ValueError("Exon start > end in %s" % (tx_id,))
                
            transcripts[tx_id]["tx_len"] += end - start
            
    return transcripts

def _tracking_parser(tracking_file, transcripts):
    with open(tracking_file, "r") as tracking_fp:
        for line_num, line in enumerate(tracking_fp):
            line_split = line.strip().split("\t")
            
            try:
                transcripts[line_split[0]]
            except KeyError:
                continue
            
            transcripts[line_split[0]]["locus_id"] = line_split[1]
            transcripts[line_split[0]]["annotation"] = line_split[2]
            transcripts[line_split[0]]["tx_classcode"] = line_split[3]
            
            for sample_data in line_split[4:]:
                if sample_data == "-": continue

                sample_split = sample_data.split("|")
                subsplit = sample_split[0].split(":")
                
                transcripts[line_split[0]]["cov"][subsplit[0]] = float(sample_split[6])

def _filter_exon_num(transcripts, min_exons):
    rejects = {}    
    
    for tx_id in transcripts.keys():
        if transcripts[tx_id]["num_exons"] < min_exons:
            rejects[tx_id] = copy.deepcopy(transcripts[tx_id])
            del transcripts[tx_id]
            
    return rejects
            
def _filter_tx_length(transcripts, min_size):
    rejects = {}    
    
    for tx_id in transcripts.keys():
        if transcripts[tx_id]["tx_len"] < min_size:
            rejects[tx_id] = copy.deepcopy(transcripts[tx_id])
            del transcripts[tx_id]
            
    return rejects
    
def _filter_coverage(transcripts, min_cov):
    rejects = {}    
    
    for tx_id in transcripts.keys():
        if max(transcripts[tx_id]["cov"].values()) < min_cov:
            rejects[tx_id] = copy.deepcopy(transcripts[tx_id])
            del transcripts[tx_id]
            
    return rejects
    
def _filter_overlaps(transcripts, filter_gtf):
    rejects = {}
    filter_tx = _gtf_parser(filter_gtf)
    
    for tx_id in transcripts.keys():
        tx_coords = zip(transcripts[tx_id]["exonStarts"].values(),
                        transcripts[tx_id]["exonEnds"].values())
        
        try:
            for filter_tx_id in filter_tx.keys():
                if filter_tx[filter_tx_id]["chrom"] != transcripts[tx_id]["chrom"]:
                    continue

                filter_coords = zip(filter_tx[filter_tx_id]["exonStarts"].values(),
                                    filter_tx[filter_tx_id]["exonEnds"].values())
                                    
                for startA, endA in tx_coords:
                    for startB, endB in filter_coords:
                        if startB <= endA <= endB or \
                           startA <= endB <= endA:
                           # OVERLAP
                           raise Exception("overlap")
        except Exception as exp_inst:
            if exp_inst.args[0] == "overlap":
                rejects[tx_id] = copy.deepcopy(transcripts[tx_id])
                del transcripts[tx_id]
            else:
                raise Exception(exp_inst.args)
            
    return rejects
    
def _filter_csf(transcripts, max_csf):
    rejects = {}    
    
    for tx_id in transcripts.keys():
        score = transcripts[tx_id]["csf"]

        if score == 'Failure("no sufficiently long ORFs found")':
            continue
        elif isinstance(score, float):
            if score > max_csf:
                rejects[tx_id] = copy.deepcopy(transcripts[tx_id])
                del transcripts[tx_id]
        else:
            raise ValueError("Unexpected value (%s) for score" % (score,))
            
    return rejects          
    
def _write_gtf(gtf, gtf_file):
    with open(gtf_file, "w") as fp:
        for tx_id in gtf.keys():
            fp.write("\n".join(gtf[tx_id]["raw"].values()))
            
def _merge_transcripts(from_tx, to_tx, note = None):
    for tx_id in from_tx.keys():
        if tx_id in to_tx:
            raise ValueError("Transcript (%s) exists in destination!" % (tx_id,))

        to_tx[tx_id] = copy.deepcopy(from_tx[tx_id])

        if note != None:
            to_tx[tx_id]["note"] = note

def _write_metadata(transcripts, outfile):
    fields = ["num_exons",
              "tx_len",
              "tx_classcode",
              "csf",
              "pfams",
              "note"]    
    
    with open(outfile, "w") as fp:
        fp.write("\t".join(["transcript_id", "max(cov)"] + fields) + "\n")
        
        for tx_id in transcripts.keys():
            fp.write("\t".join(map(str, [tx_id, max(transcripts[tx_id]["cov"].values())] + \
            [transcripts[tx_id][x] for x in fields])) + "\n")
            
def main():
    if not os.path.exists(options.output_dir):
        os.mkdir(options.output_dir)

    # parse the GTF
    transcripts = _gtf_parser(args[1])
    
    # add data from tracking file
    _tracking_parser(args[2], transcripts)
    
    all_rejects = {}
    
    rejects = _filter_exon_num(transcripts, options.min_exons)
    _write_gtf(rejects, os.path.join(options.output_dir,
                                     "rejects.too_few_exons.gtf"))
    _merge_transcripts(rejects, all_rejects, "Too few exons")
    
    rejects = _filter_tx_length(transcripts, options.min_size)
    _write_gtf(rejects, os.path.join(options.output_dir,
                                     "rejects.too_short.gtf"))
    _merge_transcripts(rejects, all_rejects, "Too short")

    rejects = _filter_coverage(transcripts, options.min_cov)
    _write_gtf(rejects, os.path.join(options.output_dir,
                                     "rejects.low_coverage.gtf"))
    _merge_transcripts(rejects, all_rejects, "Low coverage")

    for filter_num, filter_gtf in enumerate(options.filter):
        rejects = _filter_overlaps(transcripts, filter_gtf)
        _write_gtf(rejects, os.path.join(options.output_dir,
                                         "rejects.filter%s.gtf" % (filter_num,)))
        _merge_transcripts(rejects, all_rejects, "Filtered by " + filter_gtf)

    # run Pfam
    if options.exclude_pfam:
        pass
                                      
    if options.max_csf:
        # obtain CSF scores
        csf_scores = _phylocsf.score_transcripts(transcripts,
                                                 args[0],
                                                 options.num_threads)

        for tx_id, score in csf_scores.items():
            transcripts[tx_id]["csf"] = score

        # filter on CSF score
        rejects = _filter_csf(transcripts, options.max_csf)
        _write_gtf(rejects, os.path.join(options.output_dir,
                                         "rejects.high_csf.gtf"))
        _merge_transcripts(rejects, all_rejects, "CSF too high")

    _write_metadata(all_rejects, os.path.join(options.output_dir, "rejects.metadata"))
    _write_gtf(all_rejects, os.path.join(options.output_dir, "rejects.gtf"))

    _write_metadata(transcripts, os.path.join(options.output_dir, "kept.metadata"))
    _write_gtf(transcripts, os.path.join(options.output_dir, "kept.gtf"))
            
if __name__ == "__main__":
    parser = OptionParser(usage="%prog [options] <assembly> <transcripts.gtf> <transcripts.tracking>",
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
                      default=[],
                      help="remove transcripts with exons overlapping these GTFs")
                      
    # (4) Positive coding potential threshold
    parser.add_option("--max-csf",
                      dest="max_csf",
                      type="float",
                      metavar="<int>",
                      default=False,
                      help="remove transcripts with CSF scores above this value")
                      
    # (5) Known protein domain filter
    parser.add_option("--exclude-pfam",
                      dest="exclude_pfam",
                      action="store_true",
                      default=False,
                      help="remove transcripts containing Pfam domains")
                      
    options, args = parser.parse_args()

    if len(args) <> 3:
        print "Error: Incorrect number of arguments"
        parser.print_help()
        sys.exit(0)

    main()
