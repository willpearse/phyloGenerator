import pytest, sys
sys.path.append("..")
from phyloGenerator import *
Entrez.email = 'wdpearse@umn.edu'

def test_BEAST():
    #...there are not enough hours in the day to test all the combinations of BEAST calls together
    #...because then I would have checked every single way of making a phylogeny, and if you're a phylogeneticist
    #...you will (hopefully agree) just how non-helpful that would be anyway...
    seqs = list(SeqIO.parse("sequences.fasta", "fasta"))
    constraint = Phylo.read("constraint.tre", "newick")
    seqs = [[x] for x in seqs]
    muscle = alignSequences(seqs)
    
    #Default option
    assert "beast temp_BEAST.xml" == BEAST(muscle[0][0], overwrite=False, runNow=False)
    assert "beast -overwrite temp_BEAST.xml" == BEAST(muscle[0][0], runNow=False)
    current = open("temp_BEAST.xml")
    default = open("default_BEAST.xml")
    assert current.readlines() == default.readlines()
    
    #Some different settings
    assert "beast -overwrite temp_BEAST.xml" == BEAST(muscle[0][0], screenRate=100, logRate=200, chainLength=20000000, method="HKY", runNow=False)
    current = open("temp_BEAST.xml")
    hky_modified = open("hky_BEAST.xml")
    assert current.readlines() == hky_modified.readlines()

    #Constraint tree
    assert "beast -overwrite temp_BEAST.xml" == BEAST(muscle[0][0], method="HKY+GAMMA", runNow=False)
    current = open("temp_BEAST.xml")
    hky_constraint = open("hky_constraint_BEAST.xml")
    assert current.readlines() == hky_constraint.readlines()

    #Run something just to be sure...
    #...can't really check much meaningful about a Bayesian search...!
    real = BEAST(muscle[0][0], chainLength=50000)
    assert real == ('tempFinal.tre', 'temp_BEAST.xml', 'temp.trees', 'temp.log')

def test_RAxML():
    seqs = list(SeqIO.parse("sequences.fasta", "fasta"))
    constraint = Phylo.read("constraint.tre", "newick")
    seqs = [[x] for x in seqs]
    muscle = alignSequences(seqs)
    
    #Default options
    assert "raxml -s tempIn.phylip -n tempOut -f d -m GTRCAT -p " in RAxML(muscle[0][0], runNow=False)
    default = RAxML(muscle[0][0])
    assert isinstance(default, Phylo.Newick.Tree)
    assert default.is_bifurcating()

    #Restarts
    assert "raxml -s tempIn.phylip -n tempOut -f d -m GTRCAT -p " in RAxML(muscle[0][0], method='localVersion-restarts=5', runNow=False)
    assert " -N 5" in RAxML(muscle[0][0], method='localVersion-restarts=5', runNow=False)
    restarts = RAxML(muscle[0][0], method='localVersion-restarts=5')
    assert isinstance(restarts, Phylo.Newick.Tree)
    assert restarts.is_bifurcating()

    #Rapid bootstrapping
    assert "raxml -s tempIn.phylip -n tempOut -f a -m GTRCAT " in RAxML(muscle[0][0], method='localVersion-integratedBootstrap=20', runNow=False)
    assert "-N 20 -x " in RAxML(muscle[0][0], method='localVersion-integratedBootstrap=20', runNow=False)
    rapid = RAxML(muscle[0][0], method='localVersion-integratedBootstrap=20')
    assert isinstance(rapid, Phylo.Newick.Tree)
    assert rapid.is_bifurcating()

    #Constraints
    constraints = RAxML(muscle[0][0], constraint=constraint)
    assert "raxml -s tempIn.phylip -n tempOut -f d -m GTRCAT -p " in RAxML(muscle[0][0], constraint=constraint, runNow=False)
    assert " -g temp_constraint.tre" in RAxML(muscle[0][0], constraint=constraint, runNow=False)
    assert isinstance(constraints, Phylo.Newick.Tree)
    assert constraints.is_bifurcating()
    
    #Partitions
    assert "raxml -s tempIn.phylip -n tempOut -f d -m GTRCAT -p " in RAxML(muscle[0][0], method='localVersion', partitions=[0, 500, muscle[0][0].get_alignment_length()], runNow=False)
    assert " -q temp_partitions.txt" in RAxML(muscle[0][0], method='localVersion', partitions=[0, 500, muscle[0][0].get_alignment_length()], runNow=False)
    partitions = RAxML(muscle[0][0], method='localVersion', partitions=[0, 500, muscle[0][0].get_alignment_length()])
    assert isinstance(partitions, Phylo.Newick.Tree)
    assert partitions.is_bifurcating()
