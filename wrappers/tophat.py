# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrew.sczesnak@med.nyu.edu"
__date__ = "10/7/2011"
__version__ = 0.0

import os
import subprocess
import _common
import time
from lib import fastq

class TopHat():
    def __init__(self, bowtie_index, out_dir, num_threads = 1):
        self.bowtie_index = bowtie_index
        self.out_dir = out_dir
        
        self.num_threads = num_threads

    def run_multiple(self, jobs, num_procs = 1):
        pid_list = []
        
        for input_fastq, tophat_options, out_dir in jobs:
            pid = _common.fork_and_run(self.run_single,
                                       input_fastq,
                                       ([] if not tophat_options else tophat_options),
                                       (self.out_dir if not out_dir else out_dir))
            pid_list.append(pid)
            
            _common.wait_for_slot(pid_list, num_procs)
            
        _common.wait_for_slot(pid_list, 0, True)

    def run_single(self, input_fastq, tophat_options = [], out_dir = None):
        if out_dir == None:
            out_dir = self.out_dir

        if os.path.exists(out_dir):
            raise ValueError("Output directory %s already exists!" % (out_dir,))

        os.mkdir(out_dir)

        log_file = os.path.join(out_dir, "tophat.log")

        if len(input_fastq) > 2:
            raise ValueError("Maximum of two reads per sample!")
    
        _, phred_ver, readlen = fastq.validate_reads(input_fastq)
        
        if phred_ver == "phred64":
            phred_param = ["--phred64-quals"]
        elif phred_ver == "phred33":
            phred_param = []
        elif phred_ver == "solexa33":
            phred_param = ["--solexa-quals"]
        else:
            raise ValueError("Unrecognized phred version")
    
        # set segment length
        if readlen < 50:
            seg_length = int(readlen / 2)
        else:
            seg_length = 25
        
        with open(log_file, "w") as log_fp:
            tophat_cmd = [_common.TOPHAT_PATH,
                          "-p", str(self.num_threads),
                          "-o", out_dir,
                          "-z", "none",
                          "--segment-length", str(seg_length)] + \
                          phred_param + \
                          tophat_options + \
                          [self.bowtie_index] + \
                          input_fastq
                                            
            log_fp.write(" ".join(tophat_cmd) + "\n\n")
            
            tophat_proc = subprocess.Popen(tophat_cmd, stderr=log_fp)
                                            
            while tophat_proc.poll() == None:
                time.sleep(1)
            
        if os.path.exists(os.path.join(out_dir, "junctions.bed")):
            return True
        else:
            return False

def bed_to_junc(bedfile):
    juncs = []
    
    with open(bedfile, "r") as fp:
        header = fp.readline()
        
        if not header.startswith("track"):
            raise ValueError("First line in BED file not a 'track' line")

        line_num = 1
        for line in fp:
            line = line.strip()
            cols = line.split()
            line_num += 1
            if len(cols) < 12: 
                raise ValueError("Malformed line %d, missing columns" % line_num)

            chromosome = cols[0]
            orientation = cols[5]
            block_starts = [int(x) for x in cols[11].split(",")]
            block_sizes = [int(x) for x in cols[10].split(",")]
    
            left_pos = int(cols[1]) + block_starts[0] + block_sizes[0] - 1 
            right_pos = int(cols[1]) + block_starts[1] 
            juncs.append("%s\t%d\t%d\t%s" % (chromosome, left_pos, right_pos, orientation))
            
    return juncs