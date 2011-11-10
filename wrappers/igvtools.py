# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrew.sczesnak@med.nyu.edu"
__date__ = "11/2/2011"
__version__ = 0.0

import os
import subprocess
import time
import _common

class IGVTools():
    def tile(self, input_wig, output_tdf, assembly):
        if not os.path.exists(input_wig):
            raise ValueError("WIG file not found")
            
        igvtools_proc = subprocess.Popen([_common.IGVTOOLS_PATH,
                                          "tile",
                                          input_wig,
                                          output_tdf,
                                          assembly])
                                          
        while igvtools_proc.poll() == None:
            time.sleep(1)
            
        if os.path.exists(output_tdf):
            return True
        else:
            return False
