#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Unit test for pma."""

import fulenbin
import unittest
from pmanalyzer import NuclPmAnalyzer, pmstatus


class RoutineTest(unittest.TestCase):
    """RoutineTest"""

    def test_nucl_analyze_with_Y(self):
        """NuclPmAnalyzer.analyze should give the expection value"""

        pmanalyzer = NuclPmAnalyzer()
        seq = 'ATGTCGTTCTGCAGCTTCTTCGGGGGCGAGGTTTTCCAGAATCACTTTGAACCTGGCGTT'
        stdseq = 'ATGTCGTTCTGCAGCTTCTTCGGGGGCGAGGTTTTCCAGAATCACTTTGAACCTGGCGTT'
        ex_pm_status = pmstatus.Y

        pm_status, pattern = pmanalyzer.analyze(seq, stdseq)
        self.assertEqual(pm_status, ex_pm_status)
        self.assertIs(pattern, None)

    def test_nucl_analyze_with_Conserved(self):
        """NuclPmAnalyzer.analyze should give the expection value"""

        pmanalyzer = NuclPmAnalyzer()
        seq = 'ATGTCGTTCTGCAGCTTCTTCGGGGGCGAGGTTTTCCAGAATCACTTTGAACCTGGCGCT'
        stdseq = 'ATGTCGTTCTGCAGCTTCTTCGGGGGCGAGGTTTTCCAGAATCACTTTGAACCTGGCGCC'
        ex_pm_status = pmstatus.Conserved

        pm_status, pattern = pmanalyzer.analyze(seq, stdseq)
        self.assertEqual(pm_status, ex_pm_status)
        self.assertIs(pattern, None)

    def test_nucl_analyze_with_PM(self):
        """NuclPmAnalyzer.analyze should give the expection value"""

        pmanalyzer = NuclPmAnalyzer()
        seq = 'ATGTCGTTCTGCAGCTTCTTCGGGGGCGAGGTTTTCCAGAATCACTTTGAACCT'
        stdseq = 'ATGTCGTTCTGCAGCTTCTTCGGGGGCGAGGTTTTCCAGAATCACTTTGAAACT'
        ex_pm_status = pmstatus.PM

        pm_status, pattern = pmanalyzer.analyze(seq, stdseq)
        self.assertEqual(pm_status, ex_pm_status)
        self.assertIs(pattern, None)


if __name__ == '__main__':
    unittest.main()