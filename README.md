# Bio-pm
A Python3 point mutation pattern analyzing tool for nucleotide sequence.

## Installation

```
pip install biopython  #  required
python3 setup.py install
```

## Examples
### Analyze point mutation status using `pm.analyze(seq, stdseq, translate=True)`

Feed your pairwised seq and stdseq into `pm.analyze`, it will return you the correct point mutation status object.

```python
>>> import pm
>>> from pm.pattern import mutant_to_str
>>> 
>>> stdseq = "ATGGGCGCT"
>>> seq_without_pm = 'ATGGGCGCT'
>>> pm.analyze(seq_without_pm, stdseq)
<pm.status.Y object with: gaps=0, nt_pm=0, aa_pm=0, stdseq='ATGGGCGCT'>
>>> 
>>> seq_conserved = "ATGGGCGCC"
>>> pm.analyze(seq_conserved, stdseq)
<pm.status.Conserved object with: gaps=0, nt_pm=1, aa_pm=0, stdseq='ATGGGCGCC'>
>>> 
>>> seq_with_pm = 'ATGGGCGAT'
>>> pm.analyze(seq_with_pm, stdseq)
<pm.status.PM object with: gaps=0, nt_pm=1, aa_pm=1, stdseq='ATGGGCGAT'>
>>> 
>>> seq_with_gap = 'ATGGGCG-C'
>>> pm.analyze(seq_with_gap, stdseq)
<pm.status.NA object with: gaps=1, nt_pm=1, aa_pm=0, stdseq='ATGGGCGC'>
>>> 

```

### Quickly compare the point mutation status objects

Point mutation status objects have its internal order, Y > Conserved > PM > NA, when the stdseq is the same. So you can quickly do comparing between them.

```python
>>> import pm
>>> stdseq = "ATGGGCGCT"
>>> seq_without_pm = 'ATGGGCGCT'
>>> seq_conserved = "ATGGGCGCC"
>>> seq_with_pm = 'ATGGGCGAT'
>>> status_Y = pm.analyze(seq_without_pm, stdseq)
>>> status_Conserved = pm.analyze(seq_conserved, stdseq)
>>> status_PM = pm.analyze(seq_with_pm, stdseq)
>>> status_Y > status_Conserved > status_PM
True
>>> sorted([status_PM, status_Y], reverse=True)
[<pm.status.Y object with: gaps=0, nt_pm=0, aa_pm=0, stdseq='ATGGGCGCT'>, <pm.status.PM object with: gaps=0, nt_pm=1, aa_pm=1, stdseq='ATGGGCGAT'>]
>>>

```

### Generate HGVS-like mutation format basing on mutant patterns
Continues from `Quickly compare the point mutation status objects`

```python
>>> from pm.pattern import mutant_to_str
>>> status_PM
<pm.status.PM object with: gaps=0, nt_pm=1, aa_pm=1, stdseq='ATGGGCGAT'>
>>> status_PM.pattern
<pm.pattern.TranslatedPattern object at 0x2b03c9cfdc18>
>>> status_PM.pattern.list()
[((8, 'C', 'A'), (3, 'A', 'D'))]
>>> for nt_pm, aa_pm in status_PM.pattern.list():
...     print(mutant_to_str(*nt_pm) + '|' + mutant_to_str(*aa_pm))
...
8C>A|3A>D

``` 

License
----
MIT