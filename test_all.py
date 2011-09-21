# -*- coding: utf-8 -*-

import unittest
import tempfile
import filecmp

import filter_reads

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
        
    def testIlluminaVersion(self):
        # version 1.4 tests
        result = filter_reads._illumina_version(self.ver14_pe_good)
        self.assertEqual(result, filter_reads._illumina14)

        result = filter_reads._illumina_version(self.ver14_pe_bad)
        self.assertEqual(result, filter_reads._illumina14)

        result = filter_reads._illumina_version(self.ver14_pe_mixed)
        self.assertEqual(result, filter_reads._illumina14)

        result = filter_reads._illumina_version(self.ver14_pe_unordered)
        self.assertEqual(result, filter_reads._illumina14)
        
        result = filter_reads._illumina_version(self.ver14_single)
        self.assertEqual(result, filter_reads._illumina14)
        
        self.assertRaises(ValueError,
                          filter_reads._illumina_version,
                          self.ver14_bad_metadata)
            
        # version 1.8 tests
        result = filter_reads._illumina_version(self.ver18_pe_good)
        self.assertEqual(result, filter_reads._illumina18)

        result = filter_reads._illumina_version(self.ver18_pe_bad)
        self.assertEqual(result, filter_reads._illumina18)

        result = filter_reads._illumina_version(self.ver18_pe_mixed)
        self.assertEqual(result, filter_reads._illumina18)

        result = filter_reads._illumina_version(self.ver18_pe_unordered)
        self.assertEqual(result, filter_reads._illumina18)
        
        result = filter_reads._illumina_version(self.ver18_single)
        self.assertEqual(result, filter_reads._illumina18)
        
        self.assertRaises(ValueError,
                          filter_reads._illumina_version,
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
                                    
if __name__ == "__main__":
    unittest.main()