# -*- coding: utf-8 -*-
"""
Analyze PM status for nucleotide sequence.

Functions:
analyze

"""

from .pattern import parse
from .status import Y, Conserved, PM_IN_DB, Codon_optimized, PM, Null


def analyze(seq, stdseq):
    """Analyze the PM between pairwised seq and stdseq.

    Analyze the consistence of nucleotide base between a sequence and it's
    global pairwised stdandar sequence. Calculate out their PM status and 
    mutation pattern.

    Return (PM status, muation pattern).

    Args:
    seq -- nucleotide sequence
    stdseq -- glable pairwised standard sequence of seq

    """
        
    self.length = len(seq)
    if self.length != len(stdseq):
        raise ValueError("inconsistent length between seq and stdseq")
    self.seq = seq
    self.stdseq = stdseq
    self.pattern = pattern
    self.gaps = 0
    self.nt_pm = 0
    for stdv, v in pattern.mutants.values():
        if stdv + v == "--":
            raise ValueError("std_variant and variant in nucleotide "
                             "mutants cannot be gap tag '-' at the same "
                             "time")

        if stdv != '-' and v != '-':
            self.nt_pm += 1
        else:
            self.gaps += 1
    assert self.gaps == stdseq.count("-") + seq.count("-"), \
            "inconsistent gaps number between sequence and pattern."
    self.aa_pm = len([None for stdv, v in pattern.aa_mutants.values() 
                           if stdv + v != "--" and stdv != v])




class NuclPmAnalyzer(object):
    """Nucleic PM status analyzer

    """

    def __init__(self, is_pm_in_db=None, is_codon_optimized=None, with_pattern_status=None):
        """
        
        Args:
        is_pm_in_db -- A callback to check if a clone is PM_IN_DB
        is_codon_optimized -- A callback to check if a clone is 
                              Codon_optimized
        with_pattern_status -- a tuple contain pm_status, e.g. 
                               ('Conserved', 'PM_IN_DB'), indicates which 
                               status to parse its PM pattern, 
                               (default: None)
            
        """

        if with_pattern_status is None:
            with_pattern_status = ()
        self.with_pattern_status = with_pattern_status
        self.is_pm_in_db = is_pm_in_db
        self.is_codon_optimized = is_codon_optimized

    def analyze(self, seq, stdseq, gaps=None, mismatch=None):
        """Do analysis

        Return the tuple, (PM_status, PM_pattern)

        Args:
        seq -- DNA sequence
        stdseq -- standard sequence whitch sequence compring to
        gaps -- alignment gaps between seq and stdseq, if None, the amount 
                of '-' in both seq and stdseq is choosen
        mismatch -- num of mismatch nuclei bases, if None, calculated from 
                    seq and stdseq alignment

        """

        if len(stdseq) != len(seq):
            raise ValueError("length of query sequence is inconsistent with "
                            "that of subject sequence")
        if gaps is None:
            gaps = stdseq.count("-") + seq.count("-")
        else:
            if not isinstance(gaps, int):
                raise TypeError("argument gaps should be int type")
        if gaps != 0:
            return (Null, None)

        pattern_obj = mutpattern.parse(seq, stdseq, translate=True)
        if mismatch == None:
            mismatch = len([None for _, stdvar, var in pattern_obj.nt_mutants if 
                                stdvar != '-' and var != ''])
        else:
            if not isinstance(mismatch, int):
                raise TypeError("argument mismatch should be int type")

        nt_pm = mismatch
        if nt_pm == 0:  # Y
            return (Y, None)
        aa_pm = len([None for _, stdv, v in pattern_obj.aa_mutants if stdv != v])
        if aa_pm == 0:  # Conserved
            if self.is_codon_optimized is not None and \
                    self.is_codon_optimized(nt_pm, aa_pm):
                return (Codon_optimized, self._pattern(
                            Codon_optimized, pattern_obj))
            return (Conserved, self._pattern(Conserved, pattern_obj))
        else:
            if self.is_pm_in_db is not None and \
                    self.is_pm_in_db(aa_pm=aa_pm, nt_pm=nt_pm, gaps=gaps, 
                                    stdseq=stdseq, seq=seq, 
                                    pattern=pattern_obj):
                return (PM_IN_DB, self._pattern(PM_IN_DB, pattern_obj))
            return (PM, self._pattern(PM, pattern_obj))

    def _pattern(self, pm_status, pattern):
        """Return pattern depending on pm_status"""

        if pm_status not in self.with_pattern_status:
            return None
        return pattern



class ProtPmAnalyzer():
    """Proteide PM status analyzer

    Not impletement yet!

    """

    pass
