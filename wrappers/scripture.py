# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrew.sczesnak@med.nyu.edu"
__date__ = "10/7/2011"
__version__ = 0.0

import os
import subprocess
import time
import _common

class Scripture():
    def __init__(self, in_dir, out_dir, num_procs = 1):
        self.in_dir = in_dir
        self.out_dir = out_dir
        self.num_procs = num_procs
        
        self.assembly, self.chroms = _common.find_genome_files(self.in_dir)
        
        if os.path.exists(self.out_dir):
            raise ValueError("Output directory %s already exists!" % (self.out_dir,))

        os.mkdir(self.out_dir)
      
    def single_chrom(self, input_bam, chrom, scripture_options = []):
        log_file = os.path.join(self.out_dir, "scripture_" + chrom + ".log")
        out_file = os.path.join(self.out_dir, chrom + ".segments")
        size_file = os.path.join(self.in_dir, self.assembly + ".sizes")
        in_fasta = os.path.join(self.in_dir, chrom + ".fa")
        
        with open(log_file, "w") as log_fp:
            scripture_cmd = [_common.SCRIPTURE_PATH,
                             "-alignment", input_bam,
                             "-out", out_file,
                             "-sizeFile", size_file,
                             "-chr", chrom,
                             "-chrSequence", in_fasta] + \
                             scripture_options
                             
            log_fp.write(" ".join(scripture_cmd) + "\n\n")    
                             
            scripture_proc = subprocess.Popen(scripture_cmd, stderr=log_fp)
                                            
            while scripture_proc.poll() == None:
                time.sleep(1)

    def all_chrom(self, input_bam, scripture_options = []):
        out_file = os.path.join(self.out_dir, "transcripts.bed")
        pid_list = []
        
        for chrom in self.chroms:
            pid = _common.fork_and_run(self.single_chrom,
                                       input_bam,
                                       chrom,
                                       scripture_options)
                                       
            pid_list.append(pid)
            _common.wait_for_slot(pid_list, self.num_procs)
            
        _common.wait_for_slot(pid_list, 0, True)
        
        with open(out_file, "w") as final_fp:           
            for chrom in self.chroms:
                in_file = os.path.join(self.out_dir, chrom + ".segments")
                
                if not os.path.exists(in_file):
                    continue
                    
                final_fp.write(open(in_file).read())

def bed_to_gtf(in_fp, out_fp):
    for line_num, line in enumerate(in_fp):
        if not line.strip(): continue

        line_split = line.strip().split("\t")

        tx_out = [line_split[0],
                  "Scripture",
                  "transcript",
                  line_split[1],
                  line_split[2],
                  "1000",
                  ".",
                  ".",
                  "gene_id \"SCRG%06d\"; transcript_id \"SCRT%06d\";" % \
                  (line_num + 1, line_num + 1)]
                      
        out_fp.write("\t".join(map(str, tx_out)) + "\n")
                      
        for exon_num, exon_start, exon_size in zip(range(int(line_split[9])),
                                                   line_split[11].split(","),
                                                   line_split[10].split(",")):
            exon_out = [line_split[0],
                        "Scripture",
                        "exon",
                        int(line_split[1]) + int(exon_start),
                        int(line_split[1]) + int(exon_start) + int(exon_size),
                        "1000",
                        ".",
                        ".",
                        "gene_id \"SCRG%06d\"; transcript_id \"SCRT%06d\";" % \
                        (line_num + 1, line_num + 1)]
                          
            out_fp.write("\t".join(map(str, exon_out)) + "\n")