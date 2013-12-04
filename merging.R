################################################################################################
###Turning a genus-level phylogeny into a species-level phylogeny with genus-level polytomies###
################################################################################################
#Will Pearse
#9/9/2012
#This code is unchecked but *heavily* commented; I strongly advise you to check this script's output very thoroughly!

###############
###Example:####
###############

#tree <- read.tree(text='((Areplacement:10,GenusB:10):10, GenusC:10);')
#genera <- c("Areplacement","Areplacement","Areplacement","GenusB","GenusB","GenusC")
#species <- c("A robur","A ilex","A crud","B homo","B sapiens","C us")
#tree <- make.composite.with.polytomies(tree, genera, species)

###############
###FUNCTIONS###
###############
require(ape)
#Replace a given list of genera in a phylogeny with polytomes of their species
#Needs a tree, and vectors of genera and species
#NOTE: vectors should be of same length, e.g. c('Quercus', 'Quercus', 'Phocoena') and c('ilex', 'robur', 'phocoena')
make.composite.with.polytomies <- function(tree, genera, species, max.genus.age=NA){
	#Functions#
	#Binds a clade into a phylogeny, replacing a given tip, and making it ultrametric	
	#A hack to fix the ultrametric nature, but APE is a pile of shite so it's all I can do	
	bind.ultrametric.by.replacement <- function(backbone, donor, replacing.tip.label, donor.length=NA){	
		#Find the species we're replacing	
		bind.point <- which(backbone$tip.label == replacing.tip.label)	
		#Bind *badly* the two trees together	
		backbone <- bind.tree(backbone, donor, where=bind.point)	
		#Now we've bound, where's the genus located?	
		which.tip <- which(backbone$tip.label == donor$tip.label[1])	
		#What node is above the genus?	
		which.node <- backbone$edge[which(backbone$edge[,2] == which.tip),1]	
		#Where is that node in the edge.list (and thus edge.length list)	
		which.edge <- which(backbone$edge[,2] == which.node)	
		#What length is that branch?	
		tip.length <- backbone$edge.length[which.edge]
		if(is.na(donor.length)){
			#It's twice as long as it should be, so replace it
			backbone$edge.length[which.edge] <- tip.length/2
			} else {
				#Correct for the manual adjustment...
				backbone$edge.length[which.edge] <- tip.length - donor.length/2
			}
		#Return the ultrametric tree!	
		return(backbone)	
	}	
	#Make a polytomy from a given species list, optionally with a tip.length	
	make.polytomy <- function(species, tip.length=NA){	
		d.f <- data.frame(spp=factor(species))	
		polytomy <- as.phylo.formula(~spp, data=d.f)	
		if(!is.na(tip.length)) polytomy$edge.length <- rep(tip.length, length(species))	
		return(polytomy)	
	}	
	#Find the unique branch length of a tip on a given phlogeny	
	find.unique.branch.length <- function(tree, tip){	
		#What number of tip.label is the tip	
		which.tip <- which(tree$tip.label == tip)	
		#What edge is the tip.label on	
		which.edge <- which(tree$edge[,2] == which.tip)	
		#What's the edge.length for that edge?	
		tip.length <- tree$edge.length[which.edge]	
		#Return that edge.length	
		return(tip.length)	
	}	
	
	#Code#
	#The genera and species vectors *must* be characters, so convert
	genera<-as.character(genera);species<-as.character(species)
	#Go along all the genera to be replaced with polytomes
	for(genus in unique(genera)){
		#Make a vector of the species that will go into this particular polytomy
		species.to.bind <- species[genera == genus]
		#If the length of that vector is 1...
		if(length(species.to.bind) == 1){
			#...we don't need to bind anything in, so just change the tip's name
			tree$tip.label[tree$tip.label == genus] <- species.to.bind
				} else {
					#Other, find the branch length leading to the particular genus
					tip.length <- find.unique.branch.length(tree, genus)
                                        #Don't warn about edge issues unless we need to think about it...
                                        edge.warning <- NA
                                        #Set the genus to the correct age
					if(!is.na(max.genus.age)){
						if(max.genus.age*2 < tip.length){
							tip.length <- min(tip.length, max.genus.age*2)
							edge.warning <- tip.length
						}
					}
					#Make the polytomy of species, with branch lengths appropriate for that genus
					polytomy <- make.polytomy(species.to.bind, (tip.length/2))
					#Bind that polytomy in, ultrametrially, to that phylogeny
					tree <- bind.ultrametric.by.replacement(tree, polytomy, genus, edge.warning)
		}
	}
	#Return the tree
	return(tree)
}
