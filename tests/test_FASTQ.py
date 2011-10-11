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

class testFilterReads(unittest.TestCase):
    def setUp(self):
        self.ver14_pe_good = "./filter_reads/ver14_pe_good.fastq"
        self.ver14_pe_bad = "./filter_reads/ver14_pe_bad.fastq"
        self.ver14_pe_mixed = "./filter_reads/ver14_pe_mixed.fastq"
        self.ver14_pe_unordered = "./filter_reads/ver14_pe_unordered.fastq"
        self.ver14_bad_metadata = "./filter_reads/ver14_bad_metadata.fastq"
        self.ver14_single = "./filter_reads/ver14_single.fastq"

        self.ver18_pe_good = "./filter_reads/ver18_pe_good.fastq"
        self.ver18_pe_bad = "./filter_reads/ver18_pe_bad.fastq"
        self.ver18_pe_mixed = "./filter_reads/ver18_pe_mixed.fastq"
        self.ver18_pe_unordered = "./filter_reads/ver18_pe_unordered.fastq"
        self.ver18_bad_metadata = "./filter_reads/ver18_bad_metadata.fastq"
        self.ver18_single = "./filter_reads/ver18_single.fastq"
        
        self.left_tempfile = tempfile.NamedTemporaryFile("w+",delete=False)
        self.right_tempfile = tempfile.NamedTemporaryFile("w+",delete=False)
        self.orphans_tempfile = tempfile.NamedTemporaryFile("w+",delete=False)
        
    def testValidateReads(self):
        # mixed FASTQ version
        self.assertRaises(ValueError,
                          fastq.validate_reads,
                          ["./pooled_tophat/v18_100bp.fastq",
                           "./pooled_tophat/v14_100bp.fastq"])

        self.assertEqual(fastq.ILLUMINA_V14,
                         fastq.validate_reads(["./pooled_tophat/v14_36bp.fastq"]))

        self.assertEqual(fastq.ILLUMINA_V14,
                         fastq.validate_reads(["./pooled_tophat/v14_72bp.fastq"]))

        self.assertEqual(fastq.ILLUMINA_V14,
                         fastq.validate_reads(["./pooled_tophat/v14_100bp.fastq"]))

        self.assertEqual(fastq.ILLUMINA_V18,
                         fastq.validate_reads(["./pooled_tophat/v18_100bp.fastq"]))

    def testReadLen(self):
        result = fastq.fastq_readlen(self.ver14_single)
        self.assertEqual(result, 36)

        result = fastq.fastq_readlen(self.ver14_pe_mixed)
        self.assertEqual(result, 72)

        result = fastq.fastq_readlen(self.ver18_pe_mixed)
        self.assertEqual(result, 100)
        
    def testFASTQVersion(self):
        # version 1.4 tests
        result = fastq.fastq_version(self.ver14_pe_good)
        self.assertEqual(result, fastq.ILLUMINA_V14)

        result = fastq.fastq_version(self.ver14_pe_bad)
        self.assertEqual(result, fastq.ILLUMINA_V14)

        result = fastq.fastq_version(self.ver14_pe_mixed)
        self.assertEqual(result, fastq.ILLUMINA_V14)

        result = fastq.fastq_version(self.ver14_pe_unordered)
        self.assertEqual(result, fastq.ILLUMINA_V14)
        
        result = fastq.fastq_version(self.ver14_single)
        self.assertEqual(result, fastq.ILLUMINA_V14)
        
        self.assertRaises(ValueError,
                          fastq.fastq_version,
                          self.ver14_bad_metadata)
            
        # version 1.8 tests
        result = fastq.fastq_version(self.ver18_pe_good)
        self.assertEqual(result, fastq.ILLUMINA_V18)

        result = fastq.fastq_version(self.ver18_pe_bad)
        self.assertEqual(result, fastq.ILLUMINA_V18)

        result = fastq.fastq_version(self.ver18_pe_mixed)
        self.assertEqual(result, fastq.ILLUMINA_V18)

        result = fastq.fastq_version(self.ver18_pe_unordered)
        self.assertEqual(result, fastq.ILLUMINA_V18)
        
        result = fastq.fastq_version(self.ver18_single)
        self.assertEqual(result, fastq.ILLUMINA_V18)
        
        self.assertRaises(ValueError,
                          fastq.fastq_version,
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

        for line_num, line in enumerate(open(self.ver14_pe_mixed, "r")):
            if line_num % 4 == 0 and line.startswith("@"):
                line = line.strip().lstrip("@")
                mate_pair, filtered, fixed = fastq._illumina14(line)
                
                self.assertEqual(mate_pair, expected[line_num][0])
                self.assertEqual(filtered, expected[line_num][1])
                self.assertEqual(fixed, expected[line_num][2])
                
        self.assertRaises(ValueError,
                          fastq._illumina14,
                          open(self.ver14_bad_metadata, "r").readline().strip().lstrip("@"))

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

        for line_num, line in enumerate(open(self.ver18_pe_mixed, "r")):
            if line_num % 4 == 0 and line.startswith("@"):
                line = line.strip().lstrip("@")
                mate_pair, filtered, fixed = fastq._illumina18(line)
                
                self.assertEqual(mate_pair, expected[line_num][0])
                self.assertEqual(filtered, expected[line_num][1])
                self.assertEqual(fixed, expected[line_num][2])
                
        self.assertRaises(ValueError,
                          fastq._illumina18,
                          open(self.ver18_bad_metadata, "r").readline().strip().lstrip("@"))        

    def testPairedParser_v14_mixed(self):
        fastq.paired_parser(open(self.ver14_pe_mixed, "r"),
                            self.left_tempfile,
                            self.right_tempfile,
                            self.orphans_tempfile,
                            fastq._illumina14,
                            64, 25, 36)
                                    
        self.left_tempfile.seek(0)
        self.right_tempfile.seek(0)
        
        self.assertTrue(filecmp.cmp(self.left_tempfile.name,
                                    "./filter_reads/ver14_pe_mixed-LEFT.fastq"))
        self.assertTrue(filecmp.cmp(self.right_tempfile.name,
                                    "./filter_reads/ver14_pe_mixed-RIGHT.fastq"))
        
    def testPairedParser_v14_unordered(self):
        self.assertRaises(ValueError,
                          fastq.paired_parser,
                          open(self.ver14_pe_unordered, "r"),
                          self.left_tempfile,
                          self.right_tempfile,
                          self.orphans_tempfile,
                          fastq._illumina14,
                          64, 25, 36)
                          
    def testPairedParser_v18_mixed(self):
        fastq.paired_parser(open(self.ver18_pe_mixed, "r"),
                            self.left_tempfile,
                            self.right_tempfile,
                            self.orphans_tempfile,
                            fastq._illumina18,
                            33, 25, 50)

        self.left_tempfile.seek(0)
        self.right_tempfile.seek(0)
        
        self.assertTrue(filecmp.cmp(self.left_tempfile.name,
                                    "./filter_reads/ver18_pe_mixed-LEFT.fastq"))
        self.assertTrue(filecmp.cmp(self.right_tempfile.name,
                                    "./filter_reads/ver18_pe_mixed-RIGHT.fastq"))
        
    def testPairedParser_v18_unordered(self):
        self.assertRaises(ValueError,
                          fastq.paired_parser,
                          open(self.ver18_pe_unordered, "r"),
                          self.left_tempfile,
                          self.right_tempfile,
                          self.orphans_tempfile,
                          fastq._illumina18,
                          33, 25, 50)

    def testSingleParser_v14(self):
        fastq.single_parser(open(self.ver14_single, "r"),
                            self.left_tempfile,
                            fastq._illumina14,
                            64, 25, 36)

        self.left_tempfile.seek(0)

        self.assertTrue(filecmp.cmp(self.left_tempfile.name,
                                    "./filter_reads/ver14_single-FILTERED.fastq"))


    def testSingleParser_v18(self):
        fastq.single_parser(open(self.ver18_single, "r"),
                            self.left_tempfile,
                            fastq._illumina18,
                            33, 25, 18)

        self.left_tempfile.seek(0)

        self.assertTrue(filecmp.cmp(self.left_tempfile.name,
                                    "./filter_reads/ver18_single-FILTERED.fastq"))