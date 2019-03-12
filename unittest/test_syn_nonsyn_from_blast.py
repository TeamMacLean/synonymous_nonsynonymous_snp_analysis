#!/usr/bin/env python3
import sys, os
import unittest

rootPath=os.path.dirname(os.path.abspath(sys.argv[0])) + "/.."
scriptPath=rootPath + "/scripts"
sys.path.append(scriptPath)

import syn_nonsyn_from_blast

class test_syn_nonsyn_from_blast(unittest.TestCase):

    def test_set_syn_or_nonsyn(self):
        self.assertEqual(syn_nonsyn_from_blast.set_syn_or_nonsyn("HNIKK", "HNIKK"), "S")
        self.assertEqual(syn_nonsyn_from_blast.set_syn_or_nonsyn("HNIKK", "HNKKK"), "NS")

    def test_get_aminoacid_sequence(self):
        query = 'GATTGCCCATCTTCT'
        subject = 'GATTTCCCATCTTCT'
        result = [('DCPSS', 'DFPSS'), ('IAHLL', 'ISHLL'), ('LPIFF', 'FPIFF')]
        self.assertEqual(syn_nonsyn_from_blast.get_aminoacid_sequence(query, subject), result)

    def test_get_n_bases_around_mismatch(self):
        query = 'GATTGCCCATCTTCT';        subject='GATTTCCCATCTTCT';        offset = 7;        mismatches=[4];        result=[('GATTGCCCATCTTCT', 'GATTTCCCATCTTCT')]
        self.assertEqual(syn_nonsyn_from_blast.get_n_bases_around_mismatch(query, subject, offset, mismatches), result)

        query = 'AAAAAAGATTGCCCATCTTCT';        subject='TTTTTTGATTTCCCATCTTCT';        offset = 7;        mismatches=[10];        result=[('AAAGATTGCCCATC', 'TTTGATTTCCCATC')]
        self.assertEqual(syn_nonsyn_from_blast.get_n_bases_around_mismatch(query, subject, offset, mismatches), result)

    def test_get_mismatch_positions(self):

        query = 'AAATAAA'
        subject = 'AAAAAAA'
        total_mismatches = 1
        self.assertEqual(syn_nonsyn_from_blast.get_mismatch_positions(query, subject, total_mismatches), [3])



if __name__ == "__main__":
    unittest.main()
