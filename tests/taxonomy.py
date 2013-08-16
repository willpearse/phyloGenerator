#Testing sequence download functions
#Will Pearse - 2013-04-25

#Setup
import pytest, sys
sys.path.append("..")
from phyloGenerator import *
Entrez.email = 'wdpearse@umn.edu'

def test_taxonIDLookup():
    result = taxonIDLookup(1034)
    assert isinstance(result, tuple)
    assert len(result) == 4
    assert result[2] == '1034'
    assert result[3] == 'Unspecified'
    assert result[0] == 'Afipia clevelandensis'
    assert result[1][0] == "Afipia"

def test_commonLookup():
    result = commonLookup("red fox")
    assert isinstance(result, tuple)
    assert len(result) == 4
    assert result[0] == 'Vulpes vulpes'
    assert result[1][0] == "Vulpes"
    assert result[2] == '9627'
    assert result[3] == 'Vertebrate Mitochondrial'
    assert commonLookup("imaginary animal") == ()

def test_cladeSpecies():
    cladeSpec = cladeSpecies("red fox")
    commonLook = commonLookup("red fox")
    assert cladeSpec[0] == commonLook
    assert cladeSpecies("imaginary animal") == ()
    phoc = cladeSpecies("phocoena")
    assert len(phoc) > 2
    assert len(phoc[0]) == 4

def test_findLineage():
    result = findLineage("vulpes vulpes")
    assert result[0] == "vulpes vulpes"
    assert result[1] == "Vulpes"
    assert result[-1] == "cellular organisms"
    assert findLineage("imaginary animal") == ()


