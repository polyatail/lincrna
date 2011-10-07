# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrew.sczesnak@med.nyu.edu"
__date__ = "10/7/2011"
__version__ = 0.0

import os
import subprocess
import time
import _common

class Cufflinks():
    def __init__(self, out_dir, num_threads = 1):
        self.out_dir = out_dir
        self.num_threads = num_threads
        
        if os.path.exists(self.out_dir):
            raise ValueError("Output directory %s already exists!" % (self.out_dir,))

        os.mkdir(self.out_dir)        
        
    def run(self, input_bam, cufflinks_options = []):
        log_file = os.path.join(self.out_dir, "cufflinks.log")
        
        with open(log_file, "w") as log_fp:
            cufflinks_cmd = [_common.CUFFLINKS_PATH, "-q",
                             "-p", str(self.num_threads),
                             "-o", self.out_dir] + \
                             cufflinks_options + \
                             [input_bam]
        
            log_fp.write(" ".join(cufflinks_cmd) + "\n\n")                         
                             
            cufflinks_proc = subprocess.Popen(cufflinks_cmd, stderr=log_fp)
                                                
            while cufflinks_proc.poll() == None:
                time.sleep(1)