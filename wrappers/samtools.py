# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrew.sczesnak@med.nyu.edu"
__date__ = "10/7/2011"
__version__ = 0.0

import os
import subprocess
import time
import _common

class SAMTools():
    def index(self, input_bam):
        if not os.path.exists(input_bam):
            raise ValueError("BAM file not found")
        elif os.path.exists(input_bam + ".bai"):
            return True
            
        samtools_proc = subprocess.Popen([_common.SAMTOOLS_PATH,
                                          "index",
                                          input_bam])
                                          
        while samtools_proc.poll() == None:
            time.sleep(1)
            
        if os.path.exists(input_bam + ".bai"):
            return True
        else:
            return False
            
    def pileup(self, input_bam):
        if not os.path.exists(input_bam):
            raise ValueError("BAM file not found")
            
        samtools_proc = subprocess.Popen([_common.SAMTOOLS_PATH,
                                          "pileup",
                                          input_bam],
                                          stdout=subprocess.PIPE)

        samtools_reader = samtools_proc.stdout.xreadlines()

        while samtools_proc.poll() == None:
            yield samtools_reader.next()
            
    def wig(self, input_bam, output_fp, filter_less_than):
        if not os.path.exists(input_bam):
            raise ValueError("BAM file not found")
            
        samtools_proc = subprocess.Popen([_common.SAMTOOLS_PATH,
                                          "pileup",
                                          input_bam],
                                          stdout=subprocess.PIPE)

        pileup_to_wig = subprocess.Popen([_common.PILEUP_TO_WIG_PATH,
                                          str(filter_less_than)],
                                          stdin=samtools_proc.stdout,
                                          stdout=output_fp)

        while samtools_proc.poll() == None:
            time.sleep(1)
            
        while pileup_to_wig.poll() == None:
            time.sleep(1)
            
