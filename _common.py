# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrew.sczesnak@med.nyu.edu"
__date__ = "9/19/2011"
__version__ = 0.0

import os
from Bio import AlignIO

genome_mafs = {"mm9": "/comp_sync/data/foreign/ucsc/20100921_multiz30way"}
genome_fasta = {"mm9": "/comp_sync/data/foreign/ucsc/20100819_mm9_sequence"}

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