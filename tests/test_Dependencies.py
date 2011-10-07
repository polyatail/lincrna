# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrew.sczesnak@med.nyu.edu"
__date__ = "9/19/2011"
__version__ = 0.0

import sys
sys.path.append("..")
import _common
import unittest
import subprocess
import time

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

    def testScripture(self):
        scripture = subprocess.Popen([_common.SCRIPTURE_PATH],
                                     stderr=subprocess.PIPE)
        
        while scripture.poll() == None:
            time.sleep(0.1)
        
        self.assertEqual(scripture.stderr.readline().strip(), _common.SCRIPTURE_VERSION)
        
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