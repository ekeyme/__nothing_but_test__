Examples
--------

Analyze point mutation status using ``pm.analyze(seq, stdseq, translate=True)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    >>> import pm
    >>>
    >>> stdseq = 'ATGGGCGC'
    >>> seq_with_gap = 'ATGGGCG-C'
    >>> pm.analyze(seq_with_gap, stdseq)
    <pm.status.NA object with: gaps=1, nt_pm=1, aa_pm=0, stdseq='ATGGGCGC'>
    >>> 

Quickly compare between ``pm.status`` objects
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``pm.status`` objects with same stdseqs have their internal order. That is ``Y > Conserved >
PM > NA``.

.. code-block:: python

    >>> import pm
    >>>
    >>> stdseq = "ATGGGCGCT"
    >>> seq_without_pm = 'ATGGGCGCT'
    >>> seq_conserved = "ATGGGCGCC"
    >>> seq_with_pm = 'ATGGGCGAT'
    >>> status_Y = pm.analyze(seq_without_pm, stdseq)
    >>> status_Conserved = pm.analyze(seq_conserved, stdseq)
    >>> status_PM = pm.analyze(seq_with_pm, stdseq)
    >>>
    >>> status_Y > status_Conserved > status_PM
    True
    >>>

Help generate HGVS-like mutation format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*Codes continues from* `Quickly compare between pm.status objects`_

.. code-block:: python

    >>> from pm.pattern import mutant_to_str
    >>>
    >>> status_PM.pattern
    <pm.pattern.TranslatedPattern object at 0x2b03c9cfdc18>
    >>>
    >>> for nt_pm, aa_pm in status_PM.pattern.list():
    ...     print(mutant_to_str(*nt_pm) + '|' + mutant_to_str(*aa_pm))
    ...
    8C>A|3A>D