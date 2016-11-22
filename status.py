# -*- coding: utf-8 -*-
"""
PM status data structure module.

Provides a string-comparable and order-comparable data structure for 
PM status. e.g. status.Y > status.Conserved is True, also, 
status.Y > str('PM') is True.

Y
Conserved
PM_IN_DB
PM
NA

functions:
str_to_status_obj -- get singleton PM status object by conventional string PM status.

"""

SCORE_MATIX = {'Y': 8.0, 
               'Conserved': 6.0, 
               'PM_IN_DB': 4.0, 
               'PM': 2.0, 
               'NA': 0.0,
               'gaps': -16.0, 
               'aa_pm': -1.0, 
               'nt_pm': -0.2
              }
MAX_BASE = 10000000


class Y(object):
    pass


class Conserved(object):
    pass


class PM_IN_DB(object):
    pass


class PM(object):
    pass


class NA(object):
    """Base object of PM status

    """
    
    __status__ = "NA"

    def __init__(self, seq, stdseq, pattern=None, gaps=0, nt_pm=0, aa_pm=0):
        """
        
        Args:
        seq -- sequence
        stdseq -- pairwised standar sequence of seq
        pattern -- pattern object
        gaps --
        nt_pm --
        aa_pm --

        """

        self.seq = seq
        self.stdseq = stdseq
        self.length = len(seq)
        self.pattern = pattern
        self.gaps = gaps
        self.nt_pm = nt_pm
        self.aa_pm = aa_pm
        if self.length > MAX_BASE:
            raise ValueError("length must not be greater than {}".format(
                    MAX_BASE))
        self.score = self._score()

    def __str__(self):
        """Convert to string."""

        return self.__status__

    def _score(self):
        """Calculate score of this status"""

        score = SCORE_MATIX[self.__status__] \
                + SCORE_MATIX['gaps'] * self.gaps \
                + SCORE_MATIX['nt_pm'] * self.nt_pm / MAX_BASE \
                + SCORE_MATIX['aa_pm'] * self.aa_pm / MAX_BASE
        return score


def str_to_status_obj(s):
    """Get _PmStatus object by string.

    Always return the singleton PM status object. If s is not str or
    cannot find a right PM status, None will be returned.
    
    Keyword arguments:
    s -- PM status name: Y, Conserved, PM_IN_DB, Codon_optimized, PM, Null

    """

    if not isinstance(s, str):
        return None

    s = s.lower()
    if s == 'y':
        return Y
    elif s == 'conserved':
        return Conserved
    elif s == 'pm_in_db':
        return PM_IN_DB
    elif s == 'codon_optimized':
        return Codon_optimized
    elif s == 'pm':
        return PM
    elif s == 'null':
        return Null
    else:
        return None


class _PmStatus(object):
    """PM status object

    PM status base class.

    Interfaces:
    __status__ -- PM status name in object variable scope. All the subclasses should 
                impletement this attribute. 'Y', 'Conserved', 'PM_IN_DB', 
                'Codon_optimized', 'PM', 'Null'

    Methods:    
    __str__(self)
        return self.__status__

    __eq__(self, other)
        retrun self.__status__ == other.__status__

    __lt__(self, other)
        return _PmStatus.order.index(self.__status__) > _PmStatus.order.index(other.__status__)

    ...

    """

    ## Interface attribute
    # __status__ = 'Null'

    # PM status order
    order = ('Y', 'Conserved', 'PM_IN_DB', 'Codon_optimized', 'PM', 'Null')

    def __str__(self):
        """Convert to string."""

        return self.__status__

    def __eq__(self, other):
        """If pm_status is the same."""

        other = self._convert2object(other)
        if not self._is_valid_operand(other):
            return NotImplemented

        return self.__status__ == other.__status__

    def __lt__(self, other):
        """Less than."""

        other = self._convert2object(other)
        if not self._is_valid_operand(other):
            return NotImplemented

        return _PmStatus.order.index(self.__status__) > _PmStatus.order.index(other.__status__)

    def __le__(self, other):
        """Less or equal."""

        other = self._convert2object(other)
        if not self._is_valid_operand(other):
            return NotImplemented

        return _PmStatus.order.index(self.__status__) >= _PmStatus.order.index(other.__status__)

    def __gt__(self, other):
        """Great than or equal."""

        other = self._convert2object(other)
        if not self._is_valid_operand(other):
            return NotImplemented

        return _PmStatus.order.index(self.__status__) < _PmStatus.order.index(other.__status__)

    def __ge__(self, other):
        """Great or equal."""

        other = self._convert2object(other)
        if not self._is_valid_operand(other):
            return NotImplemented

        return _PmStatus.order.index(self.__status__) <= _PmStatus.order.index(other.__status__)

    def _is_valid_operand(self, other):
        """Check"""

        return isinstance(other, _PmStatus)

    def _convert2object(self, other):
        """Convert conventional string PM status into _PmStatus object.
        
        If cannot find no valid _PmStatus object, return the orignal other.

        """

        other = str_to_status_obj(other) or other
        return other


class PmStatusY(_PmStatus):
    """Y status object
    
    This object only impletement the __status__ interface attribute and inherits 
    many “rich comparison” methods from its baseclass _PmStatus.

    """
    __status__ = 'Y'


class PmStatusConserved(_PmStatus):
    """Conserved status object
    
    This object only impletement the __status__ interface attribute and inherits 
    many “rich comparison” methods from its baseclass _PmStatus.

    """
    __status__ = 'Conserved'


class PmStatusPmInDb(_PmStatus):
    """PM_IN_DB status object
    
    This object only impletement the __status__ interface attribute and inherits 
    many “rich comparison” methods from its baseclass _PmStatus.

    """
    __status__ = 'PM_IN_DB'


class PmStatusCodonOptimized(_PmStatus):
    """Codon_optimized status object
    
    This object only impletement the __status__ interface attribute and inherits 
    many “rich comparison” methods from its baseclass _PmStatus.

    """
    __status__ = 'Codon_optimized'


class PmStatusPm(_PmStatus):
    """PM PM_status object
    
    This object only impletement the __status__ interface attribute and inherits 
    many “rich comparison” methods from its baseclass _PmStatus.

    """
    __status__ = 'PM'


class PmStatusNull(_PmStatus):
    """Null PM_status object
    
    This object only impletement the __status__ interface attribute and inherits 
    many “rich comparison” methods from its baseclass _PmStatus.

    """
    __status__ = 'Null'


# PM status object
Y = PmStatusY()
Conserved =PmStatusConserved()
PM_IN_DB = PmStatusPmInDb()
Codon_optimized = PmStatusCodonOptimized()
PM = PmStatusPm()
Null = PmStatusNull()
