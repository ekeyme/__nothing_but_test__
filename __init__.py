# -*- coding: utf-8 -*-
"""
Analyze PM status for nucleotide sequence.

Functions:
analyze

"""

from .pattern import parse
from .status import Y, Conserved, PM, NA


def analyze(seq, stdseq):
    """Analyze the PM between pairwised seq and stdseq.

    Analyze the consistence of nucleotide base between a sequence and it's
    global pairwised stdandar sequence. Calculate out their PM status and 
    mutation pattern.

    Return point mutation object: Y/Convered/PM/NA.

    Args:
    seq -- nucleotide sequence
    stdseq -- glable pairwised standard sequence of seq

    """
        
    length = len(seq)
    if length != len(stdseq):
        raise ValueError("inconsistent length between seq and stdseq")
    pattern = parse(seq, stdseq, translate=True)
    gaps, nt_pm, aa_pm = 0, 0, 0
    for stdv, v in pattern.mutants.values():
        assert stdv != v, \
            "invalid pattern: consistent base exist in both stdseq and seq."
        if stdv == '-' or v == '-':
            gaps += 1
        else:
            nt_pm += 1
    assert gaps == stdseq.count("-") + seq.count("-"), \
            "inconsistent gaps number between sequence and pattern."
    aa_pm = len([None for stdv, v in pattern.aa_mutants.values() \
                    if stdv != v and stdv != '-' and v != '-'])
    if gaps > 0:
        return NA(seq, stdseq, pattern=pattern, length=length, 
                  gaps=gaps, nt_pm=nt_pm, aa_pm=aa_pm)
    if nt_pm == 0:
        return Y(seq, stdseq, pattern=pattern, length=length)
    if aa_pm == 0:
        return Conserved(seq, stdseq, pattern=pattern, length=length, nt_pm=nt_pm)
    else:
        return PM(seq, stdseq, pattern=pattern, length=length, 
                  nt_pm=nt_pm, aa_pm=aa_pm)

__all__ = (analyze, )
