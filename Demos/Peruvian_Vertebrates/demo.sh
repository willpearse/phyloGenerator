#Vertebrate example
# - data taken from a literature review by an MSc student at Imperial College London
# - thanks to Lucinda Kirkpatrick (lucinda.kirkpatrick11@imperial.ac.uk)

#DESCRIPTION
#'species.txt' - list of species
#'alignment.fasta' - COI alignment (after trimming, etc.)

#Simple search
./phyloGenerator.app/Contents/MacOS/phyloGenerator -name peruSimple -wd WHEREYOUWANTYOURFILES -email PUTINYOUROWNEMAIL -gene COI -species FILLINTHERESTOFTHISPATH/Demos/Peru_Vertebrates/species.txt -alignment mafft -phylogen beast-GTR-GAMMA
#This should look something like this, once you've filled in the bits I've put in capitals: ./phyloGenerator.app/Contents/MacOS/phyloGenerator -name peruSimple -wd /Users/bob/phyloGenTest -email bob@university.ac.uk -gene COI -species /Users/bob/phyloGenTest/Demos/Peru_Vertebrates/species.txt -alignment mafft -phylogen beast-GTR-GAMMA
#This won't turn up sequences for everything, which is fine. So, when you get some sequences, enter 'replace' mode, and type 'THOROUGH'. phyloGenerator will search for replacement sequences for you.
#You will then need to 'merge' some congeners together ('Cebus'). Enter the 'merge' mode, and do so.
#You will also need to 'trim' any long sequences (I found one ~16000bp) phyloGenerator picks out for you. Just enter 'trim' mode, specify you're dealing with mitochondrial data, and away you go. If you don't do this, you may crash your computer!
#Check the alignments for warning messages, and compare them with 'metal'.
#You can't get sequences for all of these species. So just hit enter to continue, and phyloGenerator will ignore those species. Or, try searching again with more genes (see below).

#Simple search with RAxML, then using PATHd8 to smooth
./phyloGenerator.app/Contents/MacOS/phyloGenerator -name peruSimpleSmooth -wd WHEREYOUWANTYOURFILES -email PUTINYOUROWNEMAIL -gene COI -dna FILLINTHERESTOFTHISPATH/Demos/Peru_Vertebrates/sequences.fasta -alignment mafft -phylogen raxml-integratedBootstrap=500
#Note that I've already downloaded sequences here (from the previous phyloGenerator run), and so the process is much faster
#You can't specify that you will use PATHd8 from the commandline; note you have to specify an outgroup (a species that is equally distantly related to all other species). Here there's no sensible outgroup, but it demonstrates the approach.
#If you smooth with BEAST, note that you don't need to specify an outgroup.

#Searching with more genes
./phyloGenerator.app/Contents/MacOS/phyloGenerator -name peruManygenes -wd WHEREYOUWANTYOURFILES -email PUTINYOUROWNEMAIL -gene COI,16S,cytb -species FILLINTHERESTOFTHISPATH/Demos/Peruvian_Vertebrates/species.txt -alignment mafft -phylogen raxml-integratedBootstrap=50


#In the paper, I use the BEAST tree to demonstrate how a tree can be built where, previously, one didn't exist. This is certainly not a perfect phylogeny, and two species had to be deleted due to no data (though these could be merged with other species if you wished to reduce the phylogen'ys resolution) but it can be generated in minutes, despite incomplete data on GenBank.