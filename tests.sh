#'Unit' tests
#...these should all work...

#COI butterflies everything aligned
./phyloGenerator.py -dna Demos/butterflies/sequences.fasta -name butterflyEverything -email will.pearse@gmail.com -gene COI -alignment everything -wd $HOME

#2->1 genes animals quick alignment
./phyloGenerator.py -species Demos/animals/species.txt -name animalsGeneSelect -email will.pearse@gmail.com -gene COI,ATPase -nGenes 1 -alignment quick -wd $HOME

#rbcL,matK plants quick alignment
./phyloGenerator.py -species Demos/plants/species.txt -name plantsSimpleGeneSelect -email will.pearse@gmail.com -gene rbcL,matK  -alignment quick -wd $HOME

#pre-downloaded constraint 
./phyloGenerator.py -name plantsSimpleGeneSelect -email will.pearse@gmail.com -dna Demos/plants/sequences.fasta  -alignment muscle -constr -wd $HOME


#Issues noted:
# - often the alignment statistics are the same. I've inserted a different alignment into it, run the function, and it picks up differences. The alignments quite literally are the same.