#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Unit test for pma."""

import fulenbin
import unittest
from pm import analyze
from pm.status import Y, Conserved, PM, NA
from pm.pattern import TranslatedPattern


class RoutineTest(unittest.TestCase):
    """RoutineTest"""

    def test_analyze_with_Y(self):
        """pm.analyze should give the expection value"""

        seq = 'ATGTCGTTCTGCAGCTTCTTCGGGGGCGAGGTTTTCCAGAATCACTTTGAACCTGGCGTT'
        stdseq = 'ATGTCGTTCTGCAGCTTCTTCGGGGGCGAGGTTTTCCAGAATCACTTTGAACCTGGCGTT'

        status = analyze(seq, stdseq)
        self.assertEqual(status, Y())
        self.assertEqual(status.seq, seq)
        self.assertEqual(status.stdseq, stdseq)
        self.assertEqual(status.gaps, 0)
        self.assertEqual(status.length, len(seq))
        self.assertEqual(status.nt_pm, 0)
        self.assertEqual(status.aa_pm, 0)
        self.assertIsInstance(status.pattern, TranslatedPattern)

    def test_analyze_with_Conserved(self):
        """pm.analyze should give the expection value"""

        seq = 'ATGTCGTTCTGCAGCTTCTTCGGGGGCGAGGTTTTCCAGAATCACTTTGAACCTGGCGCT'
        stdseq = 'ATGTCGTTCTGCAGCTTCTTCGGGGGCGAGGTTTTCCAGAATCACTTTGAACCTGGCGCC'

        status = analyze(seq, stdseq)
        self.assertEqual(status, Conserved(nt_pm=1))
        self.assertEqual(status.seq, seq)
        self.assertEqual(status.stdseq, stdseq)
        self.assertEqual(status.gaps, 0)
        self.assertEqual(status.length, len(seq))
        self.assertEqual(status.nt_pm, 1)
        self.assertEqual(status.aa_pm, 0)
        self.assertIsInstance(status.pattern, TranslatedPattern)

    def test_analyze_with_PM(self):
        """pm.analyze should give the expection value"""

        seq = 'ATGTCGTTCTGCAGCTTCTTCGGGGGCGAGGTTTTCCAGAATCACTTTGAACCT'
        stdseq = 'ATGTCGTTCTGCAGCTTCTTCGGGGGCGAGGTTTTCCAGAATCACTTTGAAACT'

        status = analyze(seq, stdseq)
        self.assertEqual(status, PM(nt_pm=1, aa_pm=1))
        self.assertEqual(status.seq, seq)
        self.assertEqual(status.stdseq, stdseq)
        self.assertEqual(status.gaps, 0)
        self.assertEqual(status.length, len(seq))
        self.assertEqual(status.nt_pm, 1)
        self.assertEqual(status.aa_pm, 1)
        self.assertIsInstance(status.pattern, TranslatedPattern)

    def test_analyze_with_Y(self):
        """pm.analyze should give the expection value"""

        seq = 'ATGTCGTTCTGCAGCTTCTTCGGGGGCGAGGTTTTCCAGAATCACTTTGAACCTGGC-TT'
        stdseq = 'ATGTCGTTCTGCAGCTTCTTCGGGGGCGAGGTTTTCCAGAATCACTTTGAACCTGGCGTT'

        status = analyze(seq, stdseq)
        self.assertEqual(status, NA(gaps=1))
        self.assertEqual(status.seq, seq.replace('-', ''))
        self.assertEqual(status.stdseq, stdseq)
        self.assertEqual(status.gaps, 1)
        self.assertEqual(status.length, len(seq))
        self.assertEqual(status.nt_pm, 0)
        self.assertEqual(status.aa_pm, 0)
        self.assertIsInstance(status.pattern, TranslatedPattern)


if __name__ == '__main__':
    unittest.main()