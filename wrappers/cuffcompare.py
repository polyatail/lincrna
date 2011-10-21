# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrew.sczesnak@med.nyu.edu"
__date__ = "10/7/2011"
__version__ = 0.0

import os
import subprocess
import time
import _common

class Cuffcompare():
    def __init__(self, in_dir, out_dir):
        self.in_dir = in_dir
        self.out_dir = out_dir

        if os.path.exists(self.out_dir):
            raise ValueError("Output directory %s already exists!" % (self.out_dir,))

        os.mkdir(self.out_dir)
        
    def run(self, input_gtfs, cuffcomp_options = []):
        log_file = os.path.join(self.out_dir, "cuffcompare.log")
    
        with open(log_file, "w") as log_fp:
            cuffcomp_cmd = [_common.CUFFCOMPARE_PATH, "-V",
                            "-o", self.out_dir,
                            "-s", self.in_dir] + \
                            cuffcomp_options + \
                            [" ".join(input_gtfs)]
            
            log_fp.write(" ".join(cuffcomp_cmd) + "\n\n")
            log_fp.flush()
            
            cuffcomp_proc = subprocess.Popen(cuffcomp_cmd, stderr=log_fp)
                                              
            while cuffcomp_proc.poll() == None:
                time.sleep(1)