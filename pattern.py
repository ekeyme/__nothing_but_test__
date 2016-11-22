# -*- coding: utf-8 -*-
"""
Mutation pattern analyzer for nucleotide/protein sequence

Functions:
parse(seq, stdseq, translate=False)
mutant_to_str(pos, stdvariant, variant)

Classes:
PlainPattern(mutants)
TranslatedPattern(nt_mutants, aa_mutants, assoc_dict)

Exceptions:
TranslateError

"""

from collections import OrderedDict
from Bio.Data import CodonTable
import fulenbin
from seqtool import is_stop_codon


# Translate table, used by _translate_codon
translate_table = CodonTable.ambiguous_generic_by_name["Standard"]. \
                    forward_table


def parse(seq, stdseq, translate=False):
    """Generate mutation pattern between seq and its pairwised stdseq.
    
    Return PlainPattern object if translate is False(default), else 
    return an TranslatedPattern object.
    
    Both PlainPattern and TranslatedPattern object are iterators contains
    the detail pattern dict.
    
    PlainPattern object is used for the sequence which do not need to be 
    translate into protein. PlainPattern.mutants attribute stores the detail
    pattern dict like {position: (variant_in_stdseq, variant_in_seq), ...}.

    TranslatedPattern object is used for the sequenct which need to compare 
    the variating on amino acid level. TranslatedPattern.mutants and 
    TranslatedPattern.aa_mutants store the detail nuclei and amino pattern 
    dict respectively, both in the style like 
    {position: (variant_in_stdseq, variant_in_seq), ...}.

    Args:
    seq -- sequence
    stdseq -- standard sequence
    translate -- translate model, whether to translate the DNA/RNA sequence.

    """

    seq_len = len(seq)
    if seq_len == 0:
        raise ValueError("empty sequence for seq")
    if seq_len != len(stdseq):
        raise ValueError("length must be consistent between seq and standard "
                        "seq")
    if translate and seq_len % 3 != 0:
        raise ValueError("sequence length must be triple in translate model")

    seq = seq.upper()
    stdseq = stdseq.upper()
    nt_mutants = OrderedDict()
    for i in range(0, seq_len):
        seq_base = seq[i]
        stdseq_base = stdseq[i]
        if seq_base != stdseq_base:
            nt_mutants[i+1] = (stdseq_base, seq_base)
    if not translate:
        return PlainPattern(nt_mutants)

    aa_mutants, assoc_dict = _make_translate_mutants(seq, stdseq, nt_mutants)
    return TranslatedPattern(nt_mutants, aa_mutants, assoc_dict)


def _make_translate_mutants(seq, stdseq, nt_mutant):
    """Return (aa_mutant, nt_pos2aa_assoc_dict)"""

    aa_mutant = OrderedDict()
    nt_pos2aa_assoc_dict = {}
    previous_aa_pos = -1
    for pos in nt_mutant:
        aa_pos, start, stop = _codon_slicing(pos)
        if aa_pos != previous_aa_pos:
            # translate codon in seq
            codon = seq[start:stop]
            try:
                aa = _translate_codon(codon)
            except KeyError:
                raise TranslateError("invalid codon in seq: {}".format(codon))
            # translate codon in standard sequence
            stdcodon = stdseq[start:stop]
            try:
                std_aa = _translate_codon(stdcodon)
            except KeyError:
                raise TranslateError("invalid codon in standard sequence: "
                                     "{}".format(stdcodon))

            aa_mutant[aa_pos] = (std_aa, aa)
        previous_aa_pos = aa_pos
        nt_pos2aa_assoc_dict[pos] = aa_pos
    return (aa_mutant, nt_pos2aa_assoc_dict)


def _translate_codon(codon):
    """used by _make_translate_mutant_list"""

    if is_stop_codon(codon):
        return '*'
    if '-' in codon:
        return '-'
    return translate_table[codon]


def _codon_slicing(nt_pos):
    """Return (AA_position, codon_index_start, codon_index_stop)"""

    aa_pos, reminder = divmod(nt_pos, 3)
    if reminder == 0:
        start = nt_pos - 3
    else:
        start = nt_pos - reminder
        aa_pos+=1
    return (aa_pos, start, start + 3)


def mutant_to_str(pos, stdvariant, variant, gap_string='-'):
    """Parse mutant to string

    Return HGVS-like format:

    (109, A, G) -> 109A>G
    (200, T, -) -> 200delT
    (360, -, C) -> 360insC
    (450, A, A) -> 450A=A

    Args:
    pos -- position
    stdvariant -- variant of standard sequence
    variant -- variant of align sequence
    gap_string -- a string represent a gap, default '-'

    """

    if stdvariant == variant:
        return "{}{}={}".format(pos, stdvariant, variant)

    if stdvariant == gap_string:
        mut = "{}ins{}".format(pos, variant)
    elif variant == gap_string:
        mut = "{}del{}".format(pos, stdvariant)
    else:
        mut = "{}{}>{}".format(pos, stdvariant, variant)
    return mut


class TranslatedPattern(object):
    """Translated mutation pattern data model

    Used for the sequenct which need to compare the variating on amino acid 
    .level

    Attributes:
    mutants -- original pattern dict betweenn the seq and its stdseq, e.g.
               {position: (variant_in_stdseq, variant_in_seq), ...}
    aa_mutants -- amino pattern dict betweenn the seq and its stdseq, 
                  e.g. {position: (variant_in_stdseq, variant_in_seq), ...}

    Methods:
    __init__(self, mutants, aa_mutants, assoc_dict)

    list(self):
        return [
                ((nt_pos, stdvariant, variant), (aa_pos, stdvariant, variant)), 
                ...
               ]

    """

    def __init__(self, mutants, aa_mutants, assoc_dict):
        """

        Args:
        mutants -- original pattern dict
        aa_mutants -- amino pattern dict
        assoc_dict -- 

        """

        self.mutants = mutants
        self.aa_mutants = aa_mutants
        self.assoc_dict = assoc_dict
        self._list = None

    def list(self):
        """Convert the pattern dict to flat list.

        Return [
                ((nt_pos, stdvariant, variant), (aa_pos, stdvariant, variant)), 
                ...
               ]

        """
        if self._list is None:
            self._list = []
            for nt_pos in self.mutants:
                aa_pos = self.assoc_dict[nt_pos]
                nt_item = (nt_pos, *self.mutants[nt_pos])
                aa_item = (aa_pos, *self.aa_mutants[aa_pos])
                self._list.append((nt_item, aa_item))
        return self._list


class PlainPattern(object):
    """Non translated mutation pattern data model

    Attributes:
    mutants -- a tuple contains nucl mutant list, 
               {position: (variant_standard, variant), ...}

    Methods:
    __init__(self, mutants)

    list(self):
        return [(pos, stdvariant, variant), ...]
        
    
    """

    def __init__(self, mutants):
        """
        Args:
        mutants -- mutant list

        """

        self.mutants = mutants
        self._list = None

    def list(self):
        """Convert the pattern dict to flat list.

        Return [(pos, stdvariant, variant), ...]

        """
        if self._list is None:
            self._list = []
            for nt_pos in self.mutants:
                self._list.append((nt_pos, *self.mutants[nt_pos]))
        return self._list


class TranslateError(ValueError):
    """Translation error"""

    pass


__all__ = (parse, )
