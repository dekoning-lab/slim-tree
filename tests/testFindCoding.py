import unittest
from unittest import mock
from utils import findCoding
from contextlib import redirect_stdout
import pandas as pd
import math
import os, io, pathlib

class testFindCoding(unittest.TestCase):

    def setUp(self):
        self.genome_length = 30
        self.codeFinder = findCoding.findCoding(self.genome_length)
        
    
    def test_get_coding_seqs(self):
        #Test a case with 1 coding gene - whole length - uses default from setUp
        self.assertEqual([0,self.genome_length-1], list(self.codeFinder.get_coding_regions().flatten()))
        
        #Test a case with 100% coding and 2 genes - should throw an error
        self.codeFinder.coding_ratio = 1
        self.codeFinder.gene_count = 2
        
        with self.assertRaises(SystemExit) as cm:
            with redirect_stdout(io.StringIO()) as sout:
                self.codeFinder.get_coding_seqs()
        
        self.assertEqual(cm.exception.code, 0)
        self.assertEqual(sout.getvalue(), ("Please ensure that if you have more than 1 gene, your coding ratio is not 1. Exiting.\n"))
        sout.close()
       
        
        #Test a case with 50% coding and 1 gene
        self.codeFinder.coding_ratio = 0.5
        self.codeFinder.gene_count = 1
        self.codeFinder.get_coding_seqs()
        self.assertEqual([0,(self.genome_length/2)-1], list(self.codeFinder.get_coding_regions().flatten()))
        
       
        #Test a case with 2 genes
        self.codeFinder.coding_ratio = 0.5
        self.codeFinder.gene_count = 2
        self.codeFinder.get_coding_seqs()
        self.assertEqual([0,7,22,29], list(self.codeFinder.get_coding_regions().flatten()))
        
        
        #Test a case with 3 genes
        self.codeFinder.coding_ratio = 0.5
        self.codeFinder.gene_count = 3
        self.codeFinder.get_coding_seqs()
        self.assertEqual([0,4,11,15,22,26], list(self.codeFinder.get_coding_regions().flatten()))
        
        #Test a case with a poor coding ratio
        self.codeFinder.coding_ratio = 0.75
        self.codeFinder.gene_count = 2
        self.codeFinder.get_coding_seqs()
        self.assertEqual([0,11,18,29], list(self.codeFinder.get_coding_regions().flatten()))
        
       
        #Test a case where there aren't any coding sequences
        self.codeFinder.coding_ratio = 0
        self.assertEqual(self.codeFinder.get_coding_seqs(), None)
        self.codeFinder.coding_ratio = 0.5
        self.codeFinder.gene_count = 0
        self.assertEqual(self.codeFinder.get_coding_seqs(), None)
        
        
    def test_get_coding_seqs_fasta(self):
        self.codeFinder.get_coding_seqs_fasta()
        self.assertEqual([0,self.genome_length-1], list(self.codeFinder.get_coding_regions().flatten()))
        
