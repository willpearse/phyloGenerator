#Plant example
# - data taken from a field experiment at Silwood Park
# - thanks to Kate Luckett (kathryn.luckett08@imperial.ac.uk)

#DESCRIPTION
#'species.txt' - list of plant species found at a site
#'sequences.fasta' - rbcL sequences for those species already downloaded from GenBank
#'constraint.tre' - a Newick constraint tree, generated using Phylomatic. Note that Phylomatic chops the ends off species' names (on my Mac!), so I've manually checked the species names in this file.
#'phylomatic_taxonomy.txt' - a taxonomy file to be used with Phylomatic
#'phylomatic_phylogeny.tre' - a refernce phylogeny (Davies et al. 2004) for use with Phylomatic

#Simple search
./phyloGenerator.app/Contents/MacOS/phyloGenerator -name silwoodSimple -wd . -email PUTINYOUROWNEMAIL -gene rbcL -species Demos/Silwood_Plants/species.txt -alignment mafft -phylogen beast-GTR-GAMMA

#Simple search with dated constraint tree (NOTE: using downloaded sequences now)
./phyloGenerator.app/Contents/MacOS/phyloGenerator -name silwoodConstrained -wd . -email PUTINYOUROWNEMAIL -gene rbcL -dna Demos/Silwood_Plants/sequences.fasta -alignment mafft -phylogen beast-GTR-GAMMA

#More complex search, with dated constraint tree - used in paper
./phyloGenerator.app/Contents/MacOS/phyloGenerator -name silwoodConstrained -wd . -email PUTINYOUROWNEMAIL -gene rbcL,matK -species Demos/Silwood_Plants/species.txt -alignment quick -phylogen beast-GTR-GAMMA

#Complex search but without the constraint tree - to demonstrate that the method works quite well even if you don't have data
./phyloGenerator.app/Contents/MacOS/phyloGenerator -name silwoodFree -wd . -email PUTINYOUROWNEMAIL -gene rbcL,matK -species Demos/Silwood_Plants/species.txt -alignment quick -phylogen beast-GTR-GAMMA

#...compare these with the Phylomatic constraint tree, as in the paper. Note the increased resolution in many of the clades, but that the named clades retain their dates. The ages at these nodes are drawn from a normal distribution with a mean of the stated age, but an SD of 1.
