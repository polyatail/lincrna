# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrew.sczesnak@med.nyu.edu"
__date__ = "9/19/2011"
__version__ = 0.0

import os
from Bio import AlignIO

# path to dependencies
BOWTIE_PATH = "/usr/bin/bowtie"
TOPHAT_PATH = "/usr/bin/tophat"
CUFFLINKS_PATH = "/opt/bio/OLD_VERSIONS/cufflinks-1.0.3/cufflinks"
SCRIPTURE_PATH = "/usr/bin/scripture"
CUFFCOMPARE_PATH = "/opt/bio/OLD_VERSIONS/cufflinks-1.0.3/cuffcompare"
PHYLOCSF_PATH = "/opt/bio/PhyloCSF-e346e1c/PhyloCSF"
PFAMSCAN_PATH = "/opt/bio/pfamscan/pfam_scan.pl"

# path to sequence data
genome_mafs = {"mm9": "/comp_sync/data/foreign/ucsc/20100921_multiz30way"}
genome_fasta = {"mm9": "/comp_sync/data/foreign/ucsc/20100816_mm9_sequence"}
bowtie_index = {"mm9": "/comp_sync/data/foreign/bowtie_index/mm9/basespace"}

# constants for unit tests
BOWTIE_VERSION = "bowtie version 0.12.7"
TOPHAT_VERSION = "TopHat v1.3.1"
CUFFLINKS_VERSION = "cufflinks v1.1.0"
CUFFCOMPARE_VERSION = "cuffcompare v1.1.0 (2699)"
PHYLOCSF_HELP = "usage: PhyloCSF.Linux.x86_64 parameter_set [file1 file2 ...]"
PFAMSCAN_HELP = "pfam_scan.pl: search a FASTA file against a library of Pfam HMMs"

def _load_maf(assembly):
    dirname = genome_mafs[assembly]
    found_assembly, chrom_files = _find_genome_files(dirname)

    assert assembly == found_assembly
    
    chrom_maf = {}
    
    for chrom in chrom_files:
        chrom_maf[chrom] = AlignIO.MafIO.MafIndex(os.path.join(dirname, chrom + ".mafindex"),
                                                  os.path.join(dirname, chrom + ".maf"),
                                                  assembly + "." + chrom)   
                                                  
    return chrom_maf

def _find_genome_files(dirname):
    dir_contents = os.listdir(dirname)
    
    assembly = None
    chromosome_files = []
    
    for item in dir_contents:
        if item.endswith(".sizes"):
            if assembly == None:
                assembly = item[:-6]
            else:
                raise ValueError("Multiple .sizes file in genome directory")   
        elif item.startswith("chr") and (item.endswith(".fa") or item.endswith(".maf")):
            chromosome_files.append(item.split(".")[0])
            
    if assembly == None:
        raise ValueError("Did not file any .sizes file")
            
    with open(os.path.join(dirname, assembly + ".sizes"), "r") as fp:
        all_chromosomes = [x.split("\t")[0] for x in fp]
        
    chroms_without_files = set(all_chromosomes).difference(chromosome_files)
        
    if chroms_without_files:
        raise ValueError("Chromosomes (%s) do not have sequence files" % (", ".join(chroms_without_files),))
        
    return assembly, chromosome_files