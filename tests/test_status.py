#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Unit test for pmstatus"""

import fulenbin
import unittest
from pm.status import Y, Conserved, PM, NA

class RoutineTest(unittest.TestCase):
    """Routine test."""

    def test_pm_status_gt_order(self):
        """Status should have a right order when doing gt-comparison"""

        r = Y() > Conserved() > Conserved(nt_pm=8) > PM() \
                > PM(aa_pm=5, nt_pm=10) > NA() > NA(gaps=1)
        self.assertTrue(r)

    def test_pm_status_lt_order(self):
        """Status should have a right order when doing lt-comparison"""

        r =  NA() < PM() < Conserved() < Y()
        self.assertTrue(r)

    def test_pm_status_le_order(self):
        """Status should give right value when doing le-comparison"""

        r = (Y() <= Y()) and (Conserved() <= Conserved()) and (PM() <= PM()) \
            and (NA() <= NA())
        self.assertTrue(r)

    def test_pm_status_ge_order(self):
        """Status should give right value when doing ge-comparison"""

        r = (Y() >= Y()) and (Conserved() >= Conserved()) and (PM() >= PM()) \
            and (NA() >= NA())
        self.assertTrue(r)

    def test_pm_status_eq_order(self):
        """Status should give right value when doing eq-comparison"""

        r = (Y() == Y()) and (Conserved() == Conserved()) and (PM() == PM()) \
            and (NA() == NA())
        self.assertTrue(r)

    def test_pm_status_ne_order(self):
        """Status should give right value when doing ne-comparison"""

        r = NA() != PM() != Conserved() != Y()
        self.assertTrue(r)

    def test_convert_pm_status_to_string(self):
        """Convert status object to string"""

        input_pairs = ((Y(), 'Y'), 
                       (Conserved(), 'Conserved'), 
                       (PM(), 'PM'), 
                       (NA(), 'NA'))
        for status, str_status in input_pairs:
            self.assertEqual(str(status), str_status)


if __name__ == '__main__':
    unittest.main()
