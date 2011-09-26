# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrew.sczesnak@med.nyu.edu"
__date__ = "9/19/2011"
__version__ = 0.0

import unittest
import tempfile
import filecmp
import subprocess
import time

import _common
import filter_reads
import pooled_tophat

class testDependencies(unittest.TestCase):
    def setUp(self):
        pass
    
    def testBowtie(self):
        bowtie = subprocess.Popen([_common.BOWTIE_PATH, "--version"],
                                  stdout=subprocess.PIPE)
        
        while bowtie.poll() == None:
            time.sleep(0.1)
        
        self.assertTrue(bowtie.stdout.readline().strip().endswith(_common.BOWTIE_VERSION))
        
    def testTopHat(self):
        tophat = subprocess.Popen([_common.TOPHAT_PATH, "--version"],
                                  stdout=subprocess.PIPE)
        
        while tophat.poll() == None:
            time.sleep(0.1)
        
        self.assertEqual(tophat.stdout.readline().strip(), _common.TOPHAT_VERSION)
        
    def testCufflinks(self):
        cufflinks = subprocess.Popen([_common.CUFFLINKS_PATH, "--no-update-check"],
                                     stderr=subprocess.PIPE)
        
        while cufflinks.poll() == None:
            time.sleep(0.1)
        
        self.assertEqual(cufflinks.stderr.readline().strip(), _common.CUFFLINKS_VERSION)
        
    def testCuffcompare(self):
        cuffcomp = subprocess.Popen([_common.CUFFCOMPARE_PATH],
                                     stderr=subprocess.PIPE)
        
        while cuffcomp.poll() == None:
            time.sleep(0.1)
        
        self.assertEqual(cuffcomp.stderr.readline().strip(), _common.CUFFCOMPARE_VERSION)
        
    def testPhyloCSF(self):
        phylocsf = subprocess.Popen([_common.PHYLOCSF_PATH],
                                    stdout=subprocess.PIPE)
        
        while phylocsf.poll() == None:
            time.sleep(0.1)
        
        self.assertEqual(phylocsf.stdout.readline().strip(), _common.PHYLOCSF_HELP)
        
    def testPfamScan(self):
        pfamscan = subprocess.Popen([_common.PFAMSCAN_PATH],
                                    stderr=subprocess.PIPE)
        
        while pfamscan.poll() == None:
            time.sleep(0.1)
            
        # skip first line
        _ = pfamscan.stderr.readline()

        self.assertEqual(pfamscan.stderr.readline().strip(), _common.PFAMSCAN_HELP)

class testFilterReads(unittest.TestCase):
    def setUp(self):
        self.ver14_pe_good = open("./test_data/ver14_pe_good.fastq", "r")
        self.ver14_pe_bad = open("./test_data/ver14_pe_bad.fastq", "r")
        self.ver14_pe_mixed = open("./test_data/ver14_pe_mixed.fastq", "r")
        self.ver14_pe_unordered = open("./test_data/ver14_pe_unordered.fastq", "r")
        self.ver14_bad_metadata = open("./test_data/ver14_bad_metadata.fastq", "r")
        self.ver14_single = open("./test_data/ver14_single.fastq", "r")

        self.ver18_pe_good = open("./test_data/ver18_pe_good.fastq", "r")
        self.ver18_pe_bad = open("./test_data/ver18_pe_bad.fastq", "r")
        self.ver18_pe_mixed = open("./test_data/ver18_pe_mixed.fastq", "r")
        self.ver18_pe_unordered = open("./test_data/ver18_pe_unordered.fastq", "r")
        self.ver18_bad_metadata = open("./test_data/ver18_bad_metadata.fastq", "r")
        self.ver18_single = open("./test_data/ver18_single.fastq", "r")
        
        self.left_tempfile = tempfile.NamedTemporaryFile("w+",delete=False)
        self.right_tempfile = tempfile.NamedTemporaryFile("w+",delete=False)

    def testReadLen(self):
        result = filter_reads.fastq_readlen(self.ver14_single)
        self.assertEqual(result, 36)

        result = filter_reads.readlen(self.ver14_pe_mixed)
        self.assertEqual(result, 72)

        result = filter_reads.readlen(self.ver18_pe_mixed)
        self.assertEqual(result, 100)
        
    def testFASTQVersion(self):
        # version 1.4 tests
        result = filter_reads.fastq_version(self.ver14_pe_good)
        self.assertEqual(result, filter_reads.ILLUMINA_V14)

        result = filter_reads.fastq_version(self.ver14_pe_bad)
        self.assertEqual(result, filter_reads.ILLUMINA_V14)

        result = filter_reads.fastq_version(self.ver14_pe_mixed)
        self.assertEqual(result, filter_reads.ILLUMINA_V14)

        result = filter_reads.fastq_version(self.ver14_pe_unordered)
        self.assertEqual(result, filter_reads.ILLUMINA_V14)
        
        result = filter_reads.fastq_version(self.ver14_single)
        self.assertEqual(result, filter_reads.ILLUMINA_V14)
        
        self.assertRaises(ValueError,
                          filter_reads.fastq_version,
                          self.ver14_bad_metadata)
            
        # version 1.8 tests
        result = filter_reads.fastq_version(self.ver18_pe_good)
        self.assertEqual(result, filter_reads.ILLUMINA_V18)

        result = filter_reads.fastq_version(self.ver18_pe_bad)
        self.assertEqual(result, filter_reads.ILLUMINA_V18)

        result = filter_reads.fastq_version(self.ver18_pe_mixed)
        self.assertEqual(result, filter_reads.ILLUMINA_V18)

        result = filter_reads.fastq_version(self.ver18_pe_unordered)
        self.assertEqual(result, filter_reads.ILLUMINA_V18)
        
        result = filter_reads.fastq_version(self.ver18_single)
        self.assertEqual(result, filter_reads.ILLUMINA_V18)
        
        self.assertRaises(ValueError,
                          filter_reads.fastq_version,
                          self.ver18_bad_metadata)
                          
    def testIllumina14(self):
        expected = {0:  ['1', 'N', 'HWUSI-EAS614_1:1:1:1077:5738:0'],
                    4:  ['1', 'N', 'HWUSI-EAS614_1:1:1:1077:11681:0'],
                    8:  ['1', 'N', 'HWUSI-EAS614_1:1:1:1077:8830:0'],
                    12:  ['1', 'Y', 'HWUSI-EAS614_1:1:1:998:2241:0'],
                    16:  ['1', 'Y', 'HWUSI-EAS614_1:1:1:998:7758:0'],
                    20:  ['1', 'Y', 'HWUSI-EAS614_1:1:1:998:9048:0'],
                    24:  ['2', 'N', 'HWUSI-EAS614_1:1:1:1077:5738:0'],
                    28:  ['2', 'N', 'HWUSI-EAS614_1:1:1:1077:11681:0'],
                    32:  ['2', 'N', 'HWUSI-EAS614_1:1:1:1077:8830:0'],
                    36:  ['2', 'Y', 'HWUSI-EAS614_1:1:1:998:2241:0'],
                    40:  ['2', 'Y', 'HWUSI-EAS614_1:1:1:998:7758:0'],
                    44:  ['2', 'Y', 'HWUSI-EAS614_1:1:1:998:9048:0']}

        for line_num, line in enumerate(self.ver14_pe_mixed):
            if line_num % 4 == 0 and line.startswith("@"):
                line = line.strip().lstrip("@")
                mate_pair, filtered, fixed = filter_reads._illumina14(line)
                
                self.assertEqual(mate_pair, expected[line_num][0])
                self.assertEqual(filtered, expected[line_num][1])
                self.assertEqual(fixed, expected[line_num][2])
                
        self.assertRaises(ValueError,
                          filter_reads._illumina14,
                          self.ver14_bad_metadata.readline().strip().lstrip("@"))

    def testIllumina18(self):
        expected = {0:  ['1', 'N', 'HWI-ST619:136:D06CDACXX:6:1101:1125:195244'],
                    4:  ['1', 'N', 'HWI-ST619:136:D06CDACXX:6:1101:1165:195249'],
                    8:  ['1', 'N', 'HWI-ST619:136:D06CDACXX:6:1101:1205:195249'],
                    12:  ['1', 'Y', 'HWI-ST619:136:D06CDACXX:6:1101:1125:1952'],
                    16:  ['1', 'Y', 'HWI-ST619:136:D06CDACXX:6:1101:1220:1952'],
                    20:  ['1', 'Y', 'HWI-ST619:136:D06CDACXX:6:1101:1112:1953'],
                    24:  ['2', 'N', 'HWI-ST619:136:D06CDACXX:6:1101:1125:195244'],
                    28:  ['2', 'N', 'HWI-ST619:136:D06CDACXX:6:1101:1165:195249'],
                    32:  ['2', 'N', 'HWI-ST619:136:D06CDACXX:6:1101:1205:195249'],
                    36:  ['2', 'Y', 'HWI-ST619:136:D06CDACXX:6:1101:1125:1952'],
                    40:  ['2', 'Y', 'HWI-ST619:136:D06CDACXX:6:1101:1220:1952'],
                    44:  ['2', 'Y', 'HWI-ST619:136:D06CDACXX:6:1101:1112:1953']}

        for line_num, line in enumerate(self.ver18_pe_mixed):
            if line_num % 4 == 0 and line.startswith("@"):
                line = line.strip().lstrip("@")
                mate_pair, filtered, fixed = filter_reads._illumina18(line)
                
                self.assertEqual(mate_pair, expected[line_num][0])
                self.assertEqual(filtered, expected[line_num][1])
                self.assertEqual(fixed, expected[line_num][2])
                
        self.assertRaises(ValueError,
                          filter_reads._illumina18,
                          self.ver18_bad_metadata.readline().strip().lstrip("@"))        

    def testPairedParser_v14_mixed(self):
        filter_reads._paired_parser(self.ver14_pe_mixed,
                                    self.left_tempfile,
                                    self.right_tempfile,
                                    filter_reads._illumina14)
                                    
        self.left_tempfile.seek(0)
        self.right_tempfile.seek(0)
        
        self.assertTrue(filecmp.cmp(self.left_tempfile.name,
                                    "./test_data/ver14_pe_mixed-LEFT.fastq"))
        self.assertTrue(filecmp.cmp(self.right_tempfile.name,
                                    "./test_data/ver14_pe_mixed-RIGHT.fastq"))
        
    def testPairedParser_v14_unordered(self):
        self.assertRaises(ValueError,
                          filter_reads._paired_parser,
                          self.ver14_pe_unordered,
                          self.left_tempfile,
                          self.right_tempfile,
                          filter_reads._illumina14)
                          
    def testPairedParser_v18_mixed(self):
        filter_reads._paired_parser(self.ver18_pe_mixed,
                                    self.left_tempfile,
                                    self.right_tempfile,
                                    filter_reads._illumina18)

        self.left_tempfile.seek(0)
        self.right_tempfile.seek(0)
        
        self.assertTrue(filecmp.cmp(self.left_tempfile.name,
                                    "./test_data/ver18_pe_mixed-LEFT.fastq"))
        self.assertTrue(filecmp.cmp(self.right_tempfile.name,
                                    "./test_data/ver18_pe_mixed-RIGHT.fastq"))
        
    def testPairedParser_v18_unordered(self):
        self.assertRaises(ValueError,
                          filter_reads._paired_parser,
                          self.ver18_pe_unordered,
                          self.left_tempfile,
                          self.right_tempfile,
                          filter_reads._illumina18)

    def testSingleParser_v14(self):
        filter_reads._single_parser(self.ver14_single,
                                    self.left_tempfile,
                                    filter_reads._illumina14)

        self.left_tempfile.seek(0)

        self.assertTrue(filecmp.cmp(self.left_tempfile.name,
                                    "./test_data/ver14_single-FILTERED.fastq"))


    def testSingleParser_v18(self):
        filter_reads._single_parser(self.ver18_single,
                                    self.left_tempfile,
                                    filter_reads._illumina18)

        self.left_tempfile.seek(0)

        self.assertTrue(filecmp.cmp(self.left_tempfile.name,
                                    "./test_data/ver18_single-FILTERED.fastq"))
                                    
class testPooledTopHat(unittest.TestCase):
    def setUp(self):
        pooled_tophat.parse_options(["-o", "./test_out",
                                     "-p", "30",
                                     _common.bowtie_index["mm9"],
                                     "/dev/null"])

    def testBedToJunc(self):
        expected = ['chr19\t3261580\t3262041\t-',
                    'chr19\t3265619\t3266376\t-',
                    'chr19\t3266417\t3267271\t-',
                    'chr19\t3268370\t3268721\t-', 
                    'chr19\t3268796\t3271591\t-', 
                    'chr19\t3271664\t3272923\t-', 
                    'chr19\t3272971\t3274415\t-', 
                    'chr19\t3274510\t3276821\t-', 
                    'chr19\t3276854\t3279596\t-']
                    
        result = pooled_tophat._bed_to_junc("./test_data/th_junc_good1.bed")
        
        self.assertEqual(result, expected)
        
        self.assertRaises(ValueError,
                          pooled_tophat._bed_to_junc,
                          "./test_data/th_junc_bad.bed")
    
    def testPoolJuncs(self):
        outfile = tempfile.NamedTemporaryFile("w+",delete=False)        
        result = pooled_tophat._pool_juncs(["./test_data/th_junc_good1.bed",
                                            "./test_data/th_junc_good2.bed"],
                                            outfile.name)
                                            
        outfile.seek(0)
        self.assertTrue(filecmp.cmp(outfile.name, "./test_data/th_pooled.juncs"))

    def testValidateReads(self):
        # mixed FASTQ version
        self.assertRaises(pooled_tophat._validate_reads,
                          ["./test_data/v18_100bp.fastq",
                           "./test_data/v14_100bp.fastq"])

        # mixed read lengths
        self.assertRaises(pooled_tophat._validate_reads,
                          ["./test_data/v14_36bp.fastq",
                           "./test_data/v14_72bp.fastq"])

    def testRunTH_v14_single(self):
        pooled_tophat.run_tophat("",
                                 [],
                                 "./test_data/v14_s_subset.fastq",
                                 "v14_single")

    def testRunTH_v14_paired(self):
        pooled_tophat.run_tophat("",
                                 [],
                                 "./test_data/v14_p_subset-left.fastq,./test_data/v14_p_subset-right.fastq",
                                 "v14_paired")
    
    def testRunTH_v18_single(self):
        pooled_tophat.run_tophat("",
                                 [],
                                 "./test_data/v18_s_subset.fastq",
                                 "v18_single")

    def testRunTH_v18_paired(self):
        pooled_tophat.run_tophat("",
                                 [],
                                 "./test_data/v18_p_subset-left.fastq,./test_data/v18_p_subset-right.fastq",
                                 "v18_paired")
    
    def testMain(self):
        # run all the libraries above in a pooled run
        pooled_tophat.main(["-o", "./test_out",
                            "-L", "v14s,v14p,v18s,v18p",
                            _common.bowtie_index["mm9"],
                            "./test_data/v14_s_subset.fastq",
                            "./test_data/v14_p_subset-left.fastq,./test_data/v14_p_subset-right.fastq",
                            "./test_data/v18_s_subset.fastq",
                            "./test_data/v18_p_subset-left.fastq,./test_data/v18_p_subset-right.fastq"]
    
#class testPooledTxome(unittest.TestCase):
#    def setUp(self):
#        pass
#    
#    def testRunCufflinks(self):
#        pass
#    
#    def testRunScripture(self):
#        pass
#    
#    def testRunCuffcompare(self):
#        pass
#    
#    def testIndexBAM(self):
#        pass
#    
#    def testBedToGTF(self):
#        pass
#    
#class testLincRNAClassify(unittest.TestCase):
#    def setUp(self):
#        pass
#    
#    def testGTFParser(self):
#        pass
#    
#    def testTrackingParser(self):
#        pass
#    
#    def testFilterExonNum(self):
#        pass
#    
#    def testFilterTxLen(self):
#        pass
#    
#    def testFilterCoverage(self):
#        pass
#    
#    def testFilterOverlaps(self):
#        pass
#    
#    def testWriteGTF(self):
#        pass
#    
#    def testMergeTranscripts(self):
#        pass
#    
#    def testWriteMetadata(self):
#        pass
#    
#class testPhyloCSF(unittest.TestCase):
#    def setUp(self):
#        pass
#
#    def testScoreBlock(self):
#        pass
#    
#    def testScoreTranscripts(self):
#        pass

if __name__ == "__main__":
    unittest.main()
