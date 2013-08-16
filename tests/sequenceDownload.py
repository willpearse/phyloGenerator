#Testing sequence download functions
#Will Pearse - 2013-04-25

#Setup
import pytest, sys
sys.path.append("..")
from phyloGenerator import *
Entrez.email = 'wdpearse@umn.edu'

def test_sequenceDownload():
    seq = sequenceDownload("quercus robur", ['rbcL'])
    assert isinstance(seq, tuple)
    assert seq[1] == 'rbcL'
    assert isinstance(seq[0], SeqRecord)
    seqs = sequenceDownload("quercus robur", ['rbcL'], noSeqs=2)
    assert isinstance(seq, tuple)
    assert isinstance(seqs[0], list)
    assert isinstance(seqs[0][0], SeqRecord)
    assert len(seqs[0]) == 2
    assert seqs[1] == 'rbcL'

def test_findGenes():
    x = findGenes(['quercus robur', 'quercus ilex'], [['rbcL']], download=True)
    assert x[1] == [['rbcL']]
    assert isinstance(x[0], list)
    assert isinstance(x[0][0], list)
    assert isinstance(x[0][0][0], SeqRecord)
    assert len(x[0]) == 2
    assert len(x[0][0]) == 1
    x = findGenes(['quercus robur', 'quercus ilex'], [['rbcL'], ['matK']], download=False)
    assert x[0].all()
    assert x[1] == [['rbcL'], ['matK']]
    assert isinstance(x, tuple)

