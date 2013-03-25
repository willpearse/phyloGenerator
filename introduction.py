#A tutorial-esque overview of phyloGenerator's internal functions
#Will Pearse - 25/3/2013

#This will only work if phyloGenerator.py is on your path, if not use something like os.path.append()
import phyloGenerator as pG
#You will only be able to run programs like MAFFT and RAxML if you have a 'requires' folder in your current working directory
#Run the 'setupLinux.py' script that is on the pG GitHub to generate such a folder

#Tell NCBI who you are
pG.Entrez.email='wdpearse@umn.edu'

#Load in some plant species
spp = []
with open('Demos/Silwood_Plants/species.txt') as plantSpp:
    for each in plantSpp:
        spp.append(each.strip())

#Let's make a pG 'gene object' (...rather grand name for a list really...)
#A list of lists - each list is a gene, within which are potential names for that gene
#pG uses the first name first, and calls it by that when handling it internally
genes = [['rbcL', 'ribulose bisphosphate'], ['matK']]

#Download some sequences! Note that you need to specify you want to actually download something!
#'delay' argument says how long to wait between downloads (sec), 'spacer' is how many downloads to wait until delaying
seqs = pG.findGenes(spp, genes, download=True)
#seqs[1] is the gene list you gave findGenes
#seqs[0] is a pG seqList
#seqs[0][0] is the first species' sequences
#seqs[0][0][0] is the first species' first gene (rbcL)
#seqs[0][0][1] is the first species' second gene (matK)

#Heck, let's see the lengths of these sequences (...there are more options...)
pG.checkSequenceList(seqs[0])
#Whoops! In my list, there's a really long and shitty sequence. We could find the gene in that sequence...
seqs[0][14][1] = pG.findGeneInSeq(seqs[0][14][1], ['matK'])

#Now let's align those sequences (method can be muscle, mafft, clustalo, prank, quick, everything... Read the source...)
align = pG.alignSequences(seqs[0], method='quick', nGenes=2)
#align[0] is the first gene's alignments
#align[0][0] is the first gene's first alignment (MUSCLE)
#align[0][0] is the first gene's second alignment (MAFFT)
#...

#Let's arbitrarily pick to use the MUSCLE alignments
finalAlign = [align[0][0], align[1][0]]

#Now let's make a phylogeny with these things!
pG.RAxML(finalAlign, runNow=False)
pG.BEAST(finalAlign, runNow=False)
#...something would have happened if we'd not specified runNow=False!