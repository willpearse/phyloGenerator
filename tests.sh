#'Unit' tests
#...these should all work...

#COI butterflies everything aligned
./phyloGenerator.py -dna /Users/will/Documents/code/phyloGenerator/Demos/butterflies/sequences.fasta -name butterflyEverything -email will.pearse@gmail.com -gene COI -alignment everything -wd .

#2->1 genes animals quick alignment
./phyloGenerator.py -species /Users/will/Documents/code/phyloGenerator/Demos/animals/species.txt -name animalsGeneSelect -email will.pearse@gmail.com -gene COI,ATPase -nGenes 1 -alignment quick -wd .

#rbcL,matK plants quick alignment
./phyloGenerator.py -species Demos/plants/species.txt -name plantsSimpleGeneSelect -email will.pearse@gmail.com -gene rbcL,matK  -alignment quick -wd .

#pre-downloaded constraint 
./phyloGenerator.py -name plantsSimpleGeneSelect -email will.pearse@gmail.com -dna /Users/will/Documents/code/phyloGenerator/Demos/plants/sequences.fasta  -alignment muscle -consTree newick_constraint.tre -wd .

#Use BEAST
./phyloGenerator.py -name plantsSimpleGeneSelect -email will.pearse@gmail.com -dna /Users/will/Documents/code/phyloGenerator/Demos/plants/sequences.fasta  -alignment clustalo -consTree /Users/will/Documents/code/phyloGenerator/Demos/plants/newick_constraint.tre -wd . - phylogen beast-GTR-GAMMA-chainLength=100000

#Use RAxML
./phyloGenerator.py -name plantsSimpleGeneSelect -email will.pearse@gmail.com -dna /Users/will/Documents/code/phyloGenerator/Demos/plants/sequences.fasta  -alignment clustalo -consTree /Users/will/Documents/code/phyloGenerator/Demos/plants/newick_constraint.tre -wd . -phylogen raxml-restart=2 -gene rbcL

#Use RAxML and rate smooth
./phyloGenerator.py -name plantsSimpleGeneSelect -email will.pearse@gmail.com -dna /Users/will/Documents/code/phyloGenerator/Demos/plants/sequences.fasta  -alignment clustalo -consTree -wd . -phylogen raxml-restart=2



#Issues noted:
# - often the alignment statistics are the same. I've inserted a different alignment into it, run the function, and it picks up differences. The alignments quite literally are the same.
