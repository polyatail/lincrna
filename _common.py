# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrew.sczesnak@med.nyu.edu"
__date__ = "9/19/2011"
__version__ = 0.0

import os
import sys
from Bio import AlignIO

# path to dependencies
SAMTOOLS_PATH = "/usr/bin/samtools"
BOWTIE_PATH = "/usr/bin/bowtie"
TOPHAT_PATH = "/usr/bin/tophat-andrew"
CUFFLINKS_PATH = "/opt/bio/OLD_VERSIONS/cufflinks-1.0.3/cufflinks"
SCRIPTURE_PATH = "/usr/bin/scripture"
CUFFCOMPARE_PATH = "/opt/bio/OLD_VERSIONS/cufflinks-1.0.3/cuffcompare"
PHYLOCSF_PATH = "/opt/bio/PhyloCSF-e346e1c/PhyloCSF"
PFAMSCAN_PATH = "/opt/bio/pfamscan/pfam_scan.pl"
UCSC_KENT_PATH = "/opt/bio/ucsc_kent-20100816/bin"
UCSC_HTTP_PATH = "http://hgdownload.cse.ucsc.edu/goldenPath/%(assembly)s/database/%(table)s.txt.gz"

# path to sequence data
genome_mafs = {"mm9": "/comp_sync/data/foreign/ucsc/20100921_multiz30way"}
genome_fasta = {"mm9": "/comp_sync/data/foreign/ucsc/20100816_mm9_sequence"}
bowtie_index = {"mm9": "/comp_sync/data/foreign/bowtie_index/mm9/basespace"}

# constants for unit tests
BOWTIE_VERSION = "bowtie version 0.12.7"
TOPHAT_VERSION = "TopHat v1.3.1"
CUFFLINKS_VERSION = "cufflinks v1.0.3"
SCRIPTURE_VERSION = "Using Version VPaperR3"
CUFFCOMPARE_VERSION = "cuffcompare v1.0.3 (2403)"
PHYLOCSF_HELP = "usage: PhyloCSF.Linux.x86_64 parameter_set [file1 file2 ...]"
PFAMSCAN_HELP = "pfam_scan.pl: search a FASTA file against a library of Pfam HMMs"

def _load_maf(assembly):
    dirname = genome_mafs[assembly]
    found_assembly, chrom_files = find_genome_files(dirname)

    assert assembly == found_assembly
    
    chrom_maf = {}
    
    for chrom in chrom_files:
        chrom_maf[chrom] = AlignIO.MafIO.MafIndex(os.path.join(dirname, chrom + ".mafindex"),
                                                  os.path.join(dirname, chrom + ".maf"),
                                                  assembly + "." + chrom)   
                                                  
    return chrom_maf

def find_genome_files(dirname):
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

def fork_and_run(func, *args, **kwds):
    pid = os.fork()
    
    if pid == 0:
        func(*args, **kwds)
        sys.exit(0)
    else:
        return pid
    
def wait_for_slot(pid_list, slots, finish = False):
    while len(pid_list) >= slots:
        waitpid_result = os.waitpid(0, os.WNOHANG & os.WUNTRACED)

        if waitpid_result[0] > 0:
            pid_list.remove(waitpid_result[0])
            
        if finish == True and len(pid_list) == 0:
            break