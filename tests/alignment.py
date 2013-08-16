#Testing alignment functions
#Will Pearse - 2013-04-25

#Setup
import pytest, sys
sys.path.append("..")
from phyloGenerator import *
Entrez.email = 'wdpearse@umn.edu'

def test_alignSequences():
    seqs = list(SeqIO.parse("sequences.fasta", "fasta"))
    short_seqs = seqs[0:5]
    seqs = [[x] for x in seqs]
    short_seqs = [[x] for x in short_seqs]

    #MUSCLE
    muscle = alignSequences(seqs)
    assert isinstance(muscle, list)
    assert isinstance(muscle[0], list)
    assert len(muscle) == 1
    assert len(muscle[0]) == 1
    assert isinstance(muscle[0][0], MultipleSeqAlignment)
    assert muscle[0][0].get_alignment_length() == 1671

    #MAFFT
    mafft = alignSequences(seqs, method='mafft')
    assert isinstance(mafft, list)
    assert isinstance(mafft[0], list)
    assert len(mafft) == 1
    assert len(mafft[0]) == 1
    assert isinstance(mafft[0][0], MultipleSeqAlignment)
    assert mafft[0][0].get_alignment_length() == 1671

    #Clustal-O
    clustal = alignSequences(seqs, method='clustalo')
    assert isinstance(clustal, list)
    assert isinstance(clustal[0], list)
    assert len(clustal) == 1
    assert len(clustal[0]) == 1
    assert isinstance(clustal[0][0], MultipleSeqAlignment)
    assert clustal[0][0].get_alignment_length() == 1671

    #Quick
    quick = alignSequences(seqs, method='quick')
    assert isinstance(quick, list)
    assert isinstance(quick[0], list)
    assert len(quick) == 1
    assert len(quick[0]) == 3
    assert muscle[0][0][0].seq.tostring() == quick[0][0][0].seq.tostring()
    assert mafft[0][0][0].seq.tostring() == quick[0][1][0].seq.tostring()
    assert clustal[0][0][0].seq.tostring() == quick[0][2][0].seq.tostring()

    #Prank (shorter because it takes longer)
    # - remember, prank is based on evolution, so it's not weird the alignment is longer...
    prank = alignSequences(short_seqs, method='prank')
    assert isinstance(prank, list)
    assert isinstance(prank[0], list)
    assert len(prank) == 1
    assert len(prank[0]) == 1
    assert isinstance(prank[0][0], MultipleSeqAlignment)
    assert prank[0][0].get_alignment_length() == 1679
    
    altogether = alignSequences(short_seqs, method='everything')
    assert isinstance(altogether, list)
    assert isinstance(altogether[0], list)
    assert len(altogether) == 1
    assert len(altogether[0]) == 4
    assert prank[0][0][0].seq.tostring() == altogether[0][3][0].seq.tostring()
 
