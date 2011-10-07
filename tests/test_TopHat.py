# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrew.sczesnak@med.nyu.edu"
__date__ = "9/19/2011"
__version__ = 0.0

import sys
sys.path.append("..")
import _common
import unittest
import tempfile
import filecmp
from lib import fastq

class testTopHat(unittest.TestCase):
    def setUp(self):
        pooled_tophat.parse_options(["-o", "./test_out",
                                     "-p", "30",
                                     "mm9",
                                     "/dev/null"])

    def testBedToJunc(self):
        expected = ['chr19\t3261634\t3261974\t-',
                    'chr19\t3265664\t3266335\t-',
                    'chr19\t3266458\t3267241\t-',
                    'chr19\t3268415\t3268659\t-',
                    'chr19\t3268841\t3271525\t-',
                    'chr19\t3271699\t3272887\t-',
                    'chr19\t3273034\t3274355\t-',
                    'chr19\t3274555\t3276755\t-',
                    'chr19\t3276918\t3279532\t-']                    
        result = pooled_tophat._bed_to_junc("./test_data/th_junc_good1.bed")        
        self.assertEqual(result, expected)

        expected = ['chr19\t3649599\t3652132\t-',
                    'chr19\t3649608\t3652141\t+',
                    'chr19\t3652329\t3659243\t-',
                    'chr19\t3689999\t3706645\t-',
                    'chr19\t3714624\t3715637\t+',
                    'chr19\t3715801\t3716490\t+',
                    'chr19\t3768339\t3786388\t+',
                    'chr19\t3786626\t3793067\t+',
                    'chr19\t3793214\t3796639\t+']
        result = pooled_tophat._bed_to_junc("./test_data/th_junc_good2.bed")        
        self.assertEqual(result, expected)
        
        self.assertRaises(ValueError,
                          pooled_tophat._bed_to_junc,
                          "./test_data/th_junc_bad1.bed")
        self.assertRaises(ValueError,
                          pooled_tophat._bed_to_junc,
                          "./test_data/th_junc_bad2.bed")
    
    def testPoolJuncs1(self):
        outfile = tempfile.NamedTemporaryFile("w+")        
        result = pooled_tophat._pool_juncs(["./test_data/th_junc_good1.bed",
                                            "./test_data/th_junc_good2.bed"],
                                            outfile.name)
                                            
        outfile.seek(0)
        self.assertTrue(filecmp.cmp(outfile.name, "./test_data/th_pooled1.juncs"))

    def testPoolJuncs2(self):
        outfile = tempfile.NamedTemporaryFile("w+")        
        result = pooled_tophat._pool_juncs(["./test_data/th_junc_good1.bed",
                                            "./test_data/th_junc_good1.bed"],
                                            outfile.name)
                                            
        outfile.seek(0)
        self.assertTrue(filecmp.cmp(outfile.name, "./test_data/th_pooled2.juncs"))


    
    def testRunTH_v18_single(self):
        pooled_tophat.run_tophat("",
                                 [],
                                 "./test_data/v18_s_subset.fastq",
                                 "v18_single")
                                 
        self.assertTrue(filecmp.cmp("./test_out/v18_single/accepted_hits.bam",
                                    "./test_data/pooled_tophat/v18_single/accepted_hits.bam"))
        self.assertTrue(filecmp.cmp("./test_out/v18_single/junctions.bed",
                                    "./test_data/pooled_tophat/v18_single/junctions.bed"))

    def testRunTH_v18_paired(self):
        pooled_tophat.run_tophat("",
                                 [],
                                 "./test_data/v18_p_subset-left.fastq,./test_data/v18_p_subset-right.fastq",
                                 "v18_paired")
                                 
        self.assertTrue(filecmp.cmp("./test_out/v18_paired/accepted_hits.bam",
                                    "./test_data/pooled_tophat/v18_paired/accepted_hits.bam"))
        self.assertTrue(filecmp.cmp("./test_out/v18_paired/junctions.bed",
                                    "./test_data/pooled_tophat/v18_paired/junctions.bed"))

    def testRunTH_v14_single(self):
        pooled_tophat.run_tophat("",
                                 [],
                                 "./test_data/v14_s_subset.fastq",
                                 "v14_single")
                                 
        self.assertTrue(filecmp.cmp("./test_out/v14_single/accepted_hits.bam",
                                    "./test_data/pooled_tophat/v14_single/accepted_hits.bam"))
        self.assertTrue(filecmp.cmp("./test_out/v14_single/junctions.bed",
                                    "./test_data/pooled_tophat/v14_single/junctions.bed"))

    def testRunTH_v14_paired(self):
        pooled_tophat.run_tophat("",
                                 [],
                                 "./test_data/v14_p_subset-left.fastq,./test_data/v14_p_subset-right.fastq",
                                 "v14_paired")
                                 
        self.assertTrue(filecmp.cmp("./test_out/v14_paired/accepted_hits.bam",
                                    "./test_data/pooled_tophat/v14_paired/accepted_hits.bam"))
        self.assertTrue(filecmp.cmp("./test_out/v14_paired/junctions.bed",
                                    "./test_data/pooled_tophat/v14_paired/junctions.bed"))
    
    def testFullRun(self):
        # run all the libraries above in a pooled run
        pooled_tophat.main(["-o", "./test_out",
                            "-L", "v18s,v18p,v14s,v14p",
                            "mm9",
                            "./test_data/v18_s_subset.fastq",
                            "./test_data/v18_p_subset-left.fastq,./test_data/v18_p_subset-right.fastq",
                            "./test_data/v14_s_subset.fastq",
                            "./test_data/v14_p_subset-left.fastq,./test_data/v14_p_subset-right.fastq"])
                                 
        self.assertTrue(filecmp.cmp("./test_out/prerun_v18s/accepted_hits.bam",
                                    "./test_data/pooled_tophat/prerun_v18s/accepted_hits.bam"))
        self.assertTrue(filecmp.cmp("./test_out/prerun_v18s/junctions.bed",
                                    "./test_data/pooled_tophat/prerun_v18s/junctions.bed"))
        self.assertTrue(filecmp.cmp("./test_out/prerun_v18p/accepted_hits.bam",
                                    "./test_data/pooled_tophat/prerun_v18p/accepted_hits.bam"))
        self.assertTrue(filecmp.cmp("./test_out/prerun_v18p/junctions.bed",
                                    "./test_data/pooled_tophat/prerun_v18p/junctions.bed"))
        self.assertTrue(filecmp.cmp("./test_out/prerun_v14s/accepted_hits.bam",
                                    "./test_data/pooled_tophat/prerun_v14s/accepted_hits.bam"))
        self.assertTrue(filecmp.cmp("./test_out/prerun_v14s/junctions.bed",
                                    "./test_data/pooled_tophat/prerun_v14s/junctions.bed"))
        self.assertTrue(filecmp.cmp("./test_out/prerun_v14p/accepted_hits.bam",
                                    "./test_data/pooled_tophat/prerun_v14p/accepted_hits.bam"))
        self.assertTrue(filecmp.cmp("./test_out/prerun_v14p/junctions.bed",
                                    "./test_data/pooled_tophat/prerun_v14p/junctions.bed"))

        self.assertTrue(filecmp.cmp("./test_out/v18s/accepted_hits.bam",
                                    "./test_data/pooled_tophat/v18s/accepted_hits.bam"))
        self.assertTrue(filecmp.cmp("./test_out/v18s/junctions.bed",
                                    "./test_data/pooled_tophat/v18s/junctions.bed"))
        self.assertTrue(filecmp.cmp("./test_out/v18p/accepted_hits.bam",
                                    "./test_data/pooled_tophat/v18p/accepted_hits.bam"))
        self.assertTrue(filecmp.cmp("./test_out/v18p/junctions.bed",
                                    "./test_data/pooled_tophat/v18p/junctions.bed"))
        self.assertTrue(filecmp.cmp("./test_out/v14s/accepted_hits.bam",
                                    "./test_data/pooled_tophat/v14s/accepted_hits.bam"))
        self.assertTrue(filecmp.cmp("./test_out/v14s/junctions.bed",
                                    "./test_data/pooled_tophat/v14s/junctions.bed"))
        self.assertTrue(filecmp.cmp("./test_out/v14p/accepted_hits.bam",
                                    "./test_data/pooled_tophat/v14p/accepted_hits.bam"))
        self.assertTrue(filecmp.cmp("./test_out/v14p/junctions.bed",
                                    "./test_data/pooled_tophat/v14p/junctions.bed"))