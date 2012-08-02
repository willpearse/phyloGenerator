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
./phyloGenerator.app/Contents/MacOS/phyloGenerator -name silwoodSimple -wd WHEREYOUWANTYOURFILES -email PUTINYOUROWNEMAIL -gene rbcL -species FILLINTHERESTOFTHISPATH/Demos/Silwood_Plants/species.txt -alignment mafft -phylogen beast-GTR-GAMMA
#This should look something like this, once you've filled in the bits I've put in capitals: ./phyloGenerator.app/Contents/MacOS/phyloGenerator -name silwoodSimple -wd /Users/bob/phyloGenTest -email bob@university.ac.uk -gene rbcL -species /Users/bob/phyloGenTest/Demos/Silwood_Plants/species.txt -alignment mafft -phylogen beast-GTR-GAMMA
#Consider trimming Anagallis_arvensis, and see what the effect of not setting the gene type beforehand is

#Simple search with dated constraint tree (NOTE: using downloaded sequences now)
./phyloGenerator.app/Contents/MacOS/phyloGenerator -name silwoodConstrained -wd WHEREYOUWANTYOURFILES -email PUTINYOUROWNEMAIL -gene rbcL -dna FILLINTHERESTOFTHISPATH/Demos/Silwood_Plants/sequences.fasta -alignment mafft -phylogen beast-GTR-GAMMA

#More complex search, with dated constraint tree - used in paper
./phyloGenerator.app/Contents/MacOS/phyloGenerator -name silwoodConstrainedTwoGenes -wd WHEREYOUWANTYOURFILES -email PUTINYOUROWNEMAIL -gene rbcL,matK -species FILLINTHERESTOFTHISPATH/Demos/Silwood_Plants/species.txt -alignment quick -phylogen beast-GTR-GAMMA

#Complex search but without the constraint tree - to demonstrate that the method works quite well even if you don't have data
./phyloGenerator.app/Contents/MacOS/phyloGenerator -name silwoodFree -wd WHEREYOUWANTYOURFILES -email PUTINYOUROWNEMAIL -gene rbcL,matK -species FILLINTHERESTOFTHISPATH/Demos/Silwood_Plants/species.txt -alignment quick -phylogen beast-GTR-GAMMA

#...compare these with the Phylomatic constraint tree, as in the paper. Note the increased resolution in many of the clades, but that the named clades retain their dates. The ages at these nodes are drawn from a normal distribution with a mean of the stated age, but an SD of 1.
