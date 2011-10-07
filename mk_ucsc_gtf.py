# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrew.sczesnak@med.nyu.edu"
__date__ = "10/7/2011"
__version__ = 0.0

from optparse import OptionParser
import os
import sys
from wrappers import kent

def Cufflinks_knownGene_GTF(assembly):
    kgXref = kent.fetch_ucsc_table(assembly, "kgXref")
    knownIsoforms = kent.fetch_ucsc_table(assembly, "knownIsoforms")
    knownGene = kent.fetch_ucsc_gtf(assembly, "knownGene")
    
    ucscid_to_xref = {}
    
    for line in kgXref:
        if line.startswith("#"): continue
        line_split = line.strip().split("\t")
        ucscid_to_xref[line_split[0]] = {"mRNA": line_split[1],
                                         "gene_symbol": line_split[4],
                                         "protein_id": line_split[6]}
                                         
    ucscid_to_clusterid = {}

    for line in knownIsoforms:
        if line.startswith("#"): continue
        line_split = line.strip().split("\t")
        ucscid_to_clusterid[line_split[1]] = int(line_split[0])

    output = []

    for line in knownGene:
        if not line.strip(): continue
        line_split = line.strip().split("\t", 8)
    
        ucscid = line_split[8].split("\"")[1]
        xref = ucscid_to_xref[ucscid]
        clusterid = ucscid_to_clusterid[ucscid]
    
        line = "\t".join(line_split[:8]) + "\t" + \
               " ".join(["gene_id \"CLUST%05d\";" % clusterid,
                         "transcript_id \"%s\";" % ucscid,
                         "gene_name \"%s\";" % xref["gene_symbol"],
                         "transcript_name \"%s\";" % xref["mRNA"],
                         "protein_id \"%s\";" % xref["protein_id"] if xref["protein_id"] != "" else ""])

        output.append(line)
        
    return output

def Cufflinks_mask_GTF(assembly):
    chromInfo = kent.fetch_ucsc_table(assembly, "chromInfo")

    chrM_size = 0
    masks = []

    for line in chromInfo:
        if not line.strip(): continue
        line_split = line.strip().split("\t")
        
        if line_split[0] == "chrM":
            chrM_size = int(line_split[1])
        
        masks.append(kent.fetch_ucsc_table(assembly, line_split[0] + "_rmsk"))
    
    masks.append(kent.fetch_ucsc_table(assembly, "rmsk"))
    
    rmsk = sum(masks, [])
    
    del masks
    
    output = []
    maskid = 0
    
    rmsk.append("\t".join(["."] * 5 + \
                          ["real_chrM", "0", str(chrM_size),
                           ".", "+", ".", "rRNA"]))
    rmsk.append("\t".join(["."] * 5 + \
                          ["real_chrM", "0", str(chrM_size),
                           ".", "-", ".", "rRNA"]))
    
    for line in rmsk:
        if not line.strip(): continue
        line_split = line.strip().split("\t")
        
        if line_split[11] not in ("tRNA", "rRNA"): continue
        if line_split[5] == "chrM": continue
        if line_split[5] == "real_chrM":
            line_split[5] = "chrM"
    
        maskid += 1
    
        GTF_line = [line_split[5],
                    "rmsk",
                    "exon",
                    line_split[6],
                    line_split[7],
                    "0.000000",
                    line_split[9],
                    ".",
                    "gene_id \"MASK%06d\"; transcript_id \"MASK%06d\";" % (maskid, maskid)]
                    
        output.append("\t".join(GTF_line))
        
    return output
    
def parse_options(arguments):
    global options, args

    parser = OptionParser(usage="%prog [options] <assembly>",
                          version="%prog " + str(__version__))
                          
    parser.add_option("-o",
                      dest="output_dir",
                      metavar="[./ucsc_out]",
                      default="./ucsc_out",
                      help="write output files to this directory")

    parser.add_option("--annotation",
                      dest="annotation",
                      action="store_true",
                      default=False,
                      help="generate UCSC knownGene GTF")

    parser.add_option("--mask",
                      dest="mask",
                      action="store_true",
                      default=False,
                      help="generate UCSC rRNA/tRNA mask GTF")
    
    options, args = parser.parse_args(arguments)
    
    if len(args) <> 1:
        print "Error: Incorrect number of arguments"
        parser.print_help()
        sys.exit(0)

def main(arguments=sys.argv[1:]):
    parse_options(arguments)
    
    if not os.path.exists(options.output_dir):
        os.mkdir(options.output_dir)
        
    if options.annotation:
        result = Cufflinks_knownGene_GTF(args[0])
        
        open(os.path.join(options.output_dir,
                          "knownGene.gtf"), "w").write("\n".join(result))

    if options.mask:
        result = Cufflinks_mask_GTF(args[0])
        
        open(os.path.join(options.output_dir,
                          "rRNA_tRNA_mask.gtf"), "w").write("\n".join(result))
        
if __name__ == "__main__":
    main()