#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Unit test for pmstatus"""

import fulenbin
import unittest
from pm import status 

class RoutineTest(unittest.TestCase):
    """Routine test."""

    def test_str_to_pmstatus_obj(self):
        """Function, str_to_pmstatus_obj, should return the right PM status object 
        when fed conventional string PM status"""

        pairs = ((pmstatus.Y, 'Y'), 
                    (pmstatus.Conserved, 'conserved'), (pmstatus.PM_IN_DB, 'PM_IN_DB'), 
                    (pmstatus.Codon_optimized, 'codon_optimized'), (pmstatus.PM, 'pm'), 
                    (pmstatus.Null, 'Null'))
        for status, str_status in pairs:
            self.assertIs(status, pmstatus.str_to_pmstatus_obj(str_status))

    def test_pm_status_gt_order(self):
        """PM status should have a right order when doing gt-comparison"""

        r = pmstatus.Y > pmstatus.Conserved > pmstatus.PM_IN_DB \
                > pmstatus.Codon_optimized > pmstatus.PM > pmstatus.Null
        self.assertTrue(r)

    def test_pm_status_lt_order(self):
        """PM status should have a right order when doing lt-comparison"""

        r = pmstatus.Null < pmstatus.PM < pmstatus.Codon_optimized \
                < pmstatus.PM_IN_DB < pmstatus.Conserved < pmstatus.Y 
        self.assertTrue(r)

    def test_pm_status_le_order(self):
        """PM status should give right value when doing le-comparison"""

        r = (pmstatus.Y <= pmstatus.Y) and (pmstatus.PM_IN_DB <= pmstatus.Y)
        self.assertTrue(r)

    def test_pm_status_ge_order(self):
        """PM status should give right value when doing ge-comparison"""

        r = (pmstatus.Y >= pmstatus.Y) and (pmstatus.Y >= pmstatus.PM_IN_DB)
        self.assertTrue(r)

    def test_pm_status_eq_order(self):
        """PM status should give right value when doing eq-comparison"""

        r = pmstatus.Y == pmstatus.Y 
        self.assertTrue(r)

    def test_pm_status_ne_order(self):
        """PM status should give right value when doing ne-comparison"""

        r = pmstatus.Null != pmstatus.PM != pmstatus.Codon_optimized \
                != pmstatus.PM_IN_DB != pmstatus.Conserved != pmstatus.Y 
        self.assertTrue(r)

    def test_pm_status_equal_comparing_with_str(self):
        """PM status should have the same order with conventional string PM status"""

        comparison_pairs = ((pmstatus.Y, 'Y'), 
                    (pmstatus.Conserved, 'conserved'), (pmstatus.PM_IN_DB, 'PM_IN_DB'), 
                    (pmstatus.Codon_optimized, 'codon_optimized'), (pmstatus.PM, 'pm'), 
                    (pmstatus.Null, 'Null'))
        for status, str_status in comparison_pairs:
            self.assertEqual(status, str_status)

    def test_convert_pm_status_to_string(self):
        """Convert PM status object to string"""

        input_pairs = ((pmstatus.Y, 'Y'), 
                    (pmstatus.Conserved, 'Conserved'), (pmstatus.PM_IN_DB, 'PM_IN_DB'), 
                    (pmstatus.Codon_optimized, 'Codon_optimized'), (pmstatus.PM, 'PM'), 
                    (pmstatus.Null, 'Null'))
        for status, str_status in input_pairs:
            self.assertEqual(str(status), str_status)


class ErrorTest(unittest.TestCase):
    """."""

    def test_pm_status_gt_comparing_with_invalid_operand(self):
        """pmstatus.Y should raise TypeError when cpmparing with invalid operand in gt-comparsion"""

        with self.assertRaises(TypeError):
            pmstatus.Y > 'f'

    def test_pm_status_ge_comparing_with_invalid_operand(self):
        """pmstatus.PM should raise TypeError when cpmparing with invalid operand in ge-comparsion"""

        with self.assertRaises(TypeError):
            pmstatus.PM >= 'f'

    def test_pm_status_lt_comparing_with_invalid_operand(self):
        """pmstatus.Y should raise TypeError when cpmparing with invalid operand in lt-comparsion"""

        with self.assertRaises(TypeError):
            pmstatus.Y < 'f'

    def test_pm_status_le_comparing_with_invalid_operand(self):
        """pmstatus.Y should raise TypeError when cpmparing with invalid operand in le-comparsion"""

        with self.assertRaises(TypeError):
            pmstatus.Y <= 'f'



if __name__ == '__main__':
    unittest.main()
