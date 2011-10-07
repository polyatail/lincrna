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