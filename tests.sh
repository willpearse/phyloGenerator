#'Unit' tests
#...these should all work...

./phyloGenerator.py -dna Demos/butterflies/sequences.fasta -name butterflyEverything -email will.pearse@gmail.com -gene COI -alignment everything
./phyloGenerator.py -species Demos/animals/species.txt -name animalsGeneSelect -email will.pearse@gmail.com -gene COI,ATPase -nGenes1 -alignment muscle
./phyloGenerator.py -species Demos/plants/species.txt -name plantsSimpleGeneSelect -email will.pearse@gmail.com -gene rbcL,matK -nGenes1 -alignment muscle

