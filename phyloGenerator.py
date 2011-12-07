#!/usr/bin/env python
# encoding: utf-8
"""
module.py
First attempt at putting together an automated phylogeny assembly programme for public demonstration.
This version will be funcitonal, i.e. there will be no class definitions
Created by Will Pearse on 2011-08-24.
Copyright (c) 2011 Imperial College london. All rights reserved.
"""

from Bio import Entrez #Taxonomy lookup
from Bio.Seq import Seq #Sequence manipulation
from Bio.SeqRecord import SeqRecord #Sequence manipulation
from Bio.SeqFeature import SeqFeature, FeatureLocation #Sequence manipulation
from Bio import SeqIO #Sequence manipulation
import random #Randomly select from multiple matches when downloading sequences
import numpy as np #Array and matrix sums
import subprocess, threading #Background process class
from Bio.Align.Applications import MuscleCommandline #Make a muscle commandline call
from Bio import AlignIO #Handle alignments
import os #Remove temporary files
import re #Search for files to delete
from Bio import Phylo #Load constructed phylogeny
import xml.etree.ElementTree as ET #XML parsing for BEAST
from Bio.Data import CodonTable #Codon positions for trimming sequences
import progressbar as pbar #For progress bar
import time #For waiting between sequence downloads

def taxonIDLookup(taxonID):
	#Lookup a species ID in the NCBI taxonomy database
	#INPUT: species ID
	#OUTPUT: tuple of scientific name, and (tuple) of lineage info (sorted genus-order), taxon ID, mitocondrial code of taxon
	#			-OR- empty tuple
	#TEST: commonLookup("human")
	#TEST: commonLookup("German doodlebug")
	#NOTE: order of return isn't what you'd want, probably, but done like this for historical reasons
	#NOTE: you wouldn't want to look up anything other than a species with this
	#TO-DO: Make the return type a dictionary key for ease of use
	#TO-DO: Check that the ID exists (worth it?)
	handleDownload = Entrez.efetch(db="taxonomy", id=taxonID, retmode="xml")
	resultsDownload = Entrez.read(handleDownload)
	handleDownload.close()
	scientificName = resultsDownload[0]['ScientificName']
	lineage = resultsDownload[0]['Lineage'].split("; ")
	lineage.reverse()
	lineage = tuple(lineage)
	taxId = resultsDownload[0]['TaxId']
	mitoCode = resultsDownload[0]['MitoGeneticCode']['MGCName']
	return(scientificName, lineage, taxId, mitoCode)

def commonLookup(spName):
	#Lookup a species name in the NCBI taxonomy database
	#INPUT: common species name
	#OUTPUT: tuple of scientific name, and (tuple) of lineage info (sorted genus-order), taxon ID, mitocondrial code of taxon
	#			-OR- empty tuple
	#TEST: commonLookup("human")
	#TEST: commonLookup("German doodlebug")
	#NOTE: order of return isn't what you'd want, probably, but done like this for historical reasons
	#TO-DO: Make the return type a dictionary key for ease of use
	handleSearch = Entrez.esearch(db="taxonomy", term=spName)
	resultsSearch = Entrez.read(handleSearch)
	handleSearch.close()
	if resultsSearch['IdList']:
		return taxonIDLookup(resultsSearch['IdList'][0])
	else:
		return(tuple())

def cladeSpecies(cladeName):
	#Find members of a clade in NCBI taxonomy database
	#INPUT: name of clade
	#OUTPUT: tuple of scientific name, taxon ID, 
	#TO-DO: we might not want to download all the matches ('Mammalia...), and at present this will just grab 20 by default...
	#TEST: cladeSpecies('quercus')
	#TEST: cladeSpecies('german doodlebug')
	searchTerm = cladeName + '[subtree] AND species[rank]'
	handleSearch = Entrez.esearch(db="taxonomy", term=searchTerm)
	resultsSearch = Entrez.read(handleSearch)
	if resultsSearch['IdList']:
		output = []
		for spId in resultsSearch['IdList']:
			output.append(taxonIDLookup(spId))
		return output
	else:
		return()

def findRelativeSequence(spName, geneName=None, cladeDepth=0, thorough=False, rettype='gb', titleText=None, noSeqs=1, seqChoice='random', download=True, retStart=0, retMax=20, targetLength=None, trimSeq=False, DNAtype='Standard', gapType='-', includeGenome=True, includePartial=True):
	#Download a relative's sequence (to be used if one can't be found for that target species)
	#TO-DO: Re-write this because it's horribly written in the while loop!
	genusName = spName.partition(' ')[0]
	namesTried = [genusName]
	genusDownload = sequenceDownload(spName, geneName, thorough, rettype, titleText, noSeqs, seqChoice, download, retStart, retMax, targetLength, trimSeq, DNAtype, gapType, includeGenome, includePartial)
	if genusDownload:
		return (genusDownload, namesTried)
	else:
		cladeNames = commonLookup(spName)[1]
		if cladeNames[0] != genusName:
			print "Using a new genus name:", cladeNames[0], "instead of", genusName
		else:
			del cladeNames[0]
		currentClade = 0
		while cladeDepth:
			attempt = sequenceDownload(cladeNames[currentClade], geneName, thorough, rettype, titleText, noSeqs, seqChoice, download, retStart, retMax, targetLength, trimSeq, DNAtype, gapType, includeGenome, includePartial)
			namesTried.append(cladeNames[currentClade])
			currentClade = currentClade + 1
			if attempt:
				cladeDepth = 0
			else:
				cladeDepth = cladeDepth - 1
		if attempt:
			return (attempt, namesTried)
		else:
			return (None, namesTried)

def eSearch(term, retStart=0, retMax=20, usehistory="y"):
	#Performs eSearch, by default with history on
	handle = Entrez.esearch(db="nucleotide",term=term, usehistory=usehistory, retStart=retStart, retMax=retMax)
	results = Entrez.read(handle)
	handle.close()
	return results

def eFetchSeqID(seqID, rettype='gb'):
	#Download a number of sequences given a sequenceID
	handle = Entrez.efetch(db="nucleotide", rettype=rettype, id=seqID)
	results = SeqIO.read(handle,rettype)
	handle.close()
	return results

def eFetchESearch(eSearchOutput, rettype='gb'):
	#Download a number of sequences given the output from an eSearch
	handle = Entrez.efetch(db="nucleotide", rettype=rettype, webenv=eSearchOutput['WebEnv'], query_key=eSearchOutput['QueryKey'])
	results = SeqIO.read(handle, rettype)
	handle.close()
	return results

def eSummary(seqID):
	#Download info on a particular sequence
	handle = Entrez.esummary(db="nucleotide", id=seqID)
	results = Entrez.read(handle)
	handle.close()
	return results[0]

def sequenceDownload(spName, geneName=None, thorough=False, rettype='gb', titleText=None, noSeqs=1, seqChoice='random', download=True, retStart=0, retMax=20, targetLength=None, trimSeq=False, DNAtype='Standard', gapType='-', includeGenome=True, includePartial=True):
	def dwnSeq(includeGenome, includePartial):
		#Download sequences of a given name and (optionally) gene name, to include partial matches, to include whole-genome matches, number of sequences, how to choose among multiple hits
		#INPUT: scientific name, optionally gene name, number of sequences, exclude partial matches, exclude whole-genome matches, multiple hits choice method, whether to actually download sequences
		#OUTPUT: Seq object(s) of that DNA region
		#			-OR- empty tuple
		#			-OR- tuple of WebEnv, QueryKey and no. hits if not downloading
		#TEST: sequenceDownload("Homo sapiens")
		#TEST: sequenceDownload("Homo sapiens", includePartial=False)
		#TEST: sequenceDownload("Homo sapiens", noSeqs=2)
		#TEST: sequenceDownload("Homo sapiens", includePartial=False, geneName="COI")
		#TEST: sequenceDownload("Homo sapiens", includeGenome=False, geneName="COI")
		#TEST: sequenceDownload("Homo sapiens", download=False)
		#TEST: sequenceDownload("German doodlebug")
		#TO-DO: include more ways to deal with multiple hits
		#TO-DO: allow custom searches to be specified
		#TO-DO: neaten the length stuff with a wrapper function - you're doing the same thing five times!
		searchTerm = spName + "[Organism]"
		if geneName: searchTerm = "(" + searchTerm + ") AND " + geneName + "[Gene Name]"
		if titleText: searchTerm = "(" + searchTerm + ") AND " + titleText + " [Title]"
		if not includePartial: searchTerm = "(" + searchTerm + ") NOT " + "partial [Title]"
		if not includeGenome: searchTerm = "(" + searchTerm + ") NOT " + "genome [Title]"
		firstSearch = eSearch(searchTerm, retStart, retMax)
		if firstSearch:
			if int(firstSearch['Count']) == noSeqs:
				return eFetchESearch(firstSearch)
			elif int(firstSearch['Count']) > noSeqs:
				if seqChoice == 'random':
					chosenIDs = random.sample(firstSearch['IdList'], noSeqs)
					if len(chosenIDs) > 1:
						return [eFetchSeqID(x) for x in chosenIDs]
					else:
						return eFetchSeqID(chosenIDs, rettype=rettype)
				elif seqChoice == 'minLength':
					if noSeqs > 1: raise RuntimeError("You can't return more than one median-length sequence...")
					currentMinLength = 1000000
					currentBest = 0
					for index, seq in enumerate(firstSearch['IdList']):
						currentLength = eSummary(seq)['Length']
						if currentLength < currentMinLength:
							currentMinLength = index
							currentBest = 0
					return eFetchSeqID(firstSearch['IdList'][currentBest], rettype=rettype)
				elif seqChoice == 'targetLength' and targetLength:
					if noSeqs > 1: raise RuntimeError("You can't return more than one best-length sequence...")
					currentMinLength = 1000000
					currentBest = 0
					for index, seq in enumerate(firstSearch['IdList']):
						currentLength = abs(targetLength - eSummary(seq)['Length'])
						if currentLength < currentMinLength:
							currentMinLength = currentLength
							currentBest = 0
					return eFetchSeqID(firstSearch['IdList'][currentBest], rettype=rettype)
				elif seqChoice == 'medianLength':
					if noSeqs > 1: raise RuntimeError("You can't return more than one best-length sequence...")
					lengths = []
					for seq in firstSearch['IdList']:
						lengths.append(eSummary(seq)['Length'])
					medianLength = np.median(lengths)
					#This next bit is necessary because you can't search for a median = a median can not be in the sample... (e.g. it could be floating point)
					currentMinLength = 1000000
					currentBest = 0
					for index, length in enumerate(lengths):
						if abs(length - medianLength) < currentMinLength:
							currentBest = firstSearch['IdList'][index]
					return eFetchSeqID(currentBest, rettype=rettype)
				else:
					raise RuntimeError("Unrecognised sequence selection method")
			else:
				return ()
		else:
			if int(firstSearch['Count']):
				return (firstSearch['WebEnv'], firstSearch['QueryKey'], int(firstSearch['Count']))
			else:
				return ()
	
	if thorough:
		seq = dwnSeq(includeGenome=False, includePartial=False)
		if seq:
			return seq
		else:
			seq = dwnSeq(includeGenome=True, includePartial=False)
			if seq:
				seq = findGeneInSeq(seq, geneName, trimSeq=trimSeq, DNAtype=DNAtype, gapType=gapType)
				return seq
			else:
				seq = dwnSeq(includeGenome=True, includePartial=True)
				if seq:
					return seq
				else:
					return ()
	else:
		return dwnSeq(includeGenome, includePartial)

def findGenes(speciesList, geneNames, download=False, titleText=None, targetNoGenes=None, noSeqs=1, includePartial=True, includeGenome=True, seqChoice='random', verbose=True, thorough=False):
	#Given a set of species, find the best genes to use to yield complete coverage of the group
	#Allows you to select the number of genes you want, or just checks all of them for you
	#OUTPUT: NumPy array of the species/genes matches, in the order each was given
	#		-OR- the genes you should use
	#ASSUMES that all these species have been checked for validity beforehand (i.e. that they're worth searching for)
	#TO-DO: Do you need to implement some kind of WebEnv/QueryKey method for this?
	#TO-DO: Let the user select the top genes after they've done a general search
	#TO-DO: Some kind of progress bar
	#TEST: findGenes(['Quercus robur', 'Quercus ilex', 'Pinus sylvaticus'], ['rbcL', 'ITS1', 'matK', 'ITS2'])
	#TEST: findGenes(['Quercus robur', 'Quercus ilex', 'Pinus sylvaticus'], ['rbcL', 'ITS1', 'matK', 'ITS2'], targetNoGenes=2)
	#TEST: findGenes(['Quercus robur', 'Quercus ilex', 'Pinus sylvaticus'], ['rbcL', 'ITS1'], targetNoGenes=2)
	if type(geneNames) is list:
		if targetNoGenes is len(geneNames):
			raise RuntimeError("Number of genes and target number of genes are the same.")
		def findBestGene(foundSeqsArray):
			geneHits = foundSeqsArray.sum(axis=0)
			for i in range(len(geneHits)):
				if geneHits[i] == max(geneHits):
					return i
		
		#Download number of genes and histories for each species
		searchResults = []
		if download:
			foundSeqs = []
		else:
			foundSeqs = np.zeros((len(speciesList), len(geneNames)), int)
		if verbose:	progBar = pbar.ProgressBar(widgets=[pbar.Percentage(), pbar.Bar()], maxval=len(speciesList)).start()
		for i in range(len(speciesList)):
			if verbose: progBar.update(i+1)
			speciesGenes = []
			for k in range(len(geneNames)):
				sequence = sequenceDownload(speciesList[i], geneNames[k], titleText=titleText, noSeqs=noSeqs, includePartial=includePartial, includeGenome=includeGenome, seqChoice=seqChoice, download=download, thorough=thorough)
				if download:
					speciesGenes.append(sequence)
				elif sequence:
					foundSeqs[i,k] = 1
			if download:
				foundSeqs.append(speciesGenes)
		if verbose: progBar.finish()
		if targetNoGenes:
			currentFoundSeqs = foundSeqs
			currentGeneNames = geneNames
			bestGenes = []
			for each in range(targetNoGenes):
				bestGeneIndex = findBestGene(currentFoundSeqs)
				bestGenes.append(currentGeneNames[bestGeneIndex])
				currentFoundSeqs = np.delete(currentFoundSeqs, bestGeneIndex, axis=1)
				del currentGeneNames[bestGeneIndex]
			return bestGenes
		else:
			return foundSeqs
	else:
		output = []
		for species in speciesList:
			sequence = sequenceDownload(species, geneNames, titleText=titleText, noSeqs=noSeqs, includePartial=includePartial, includeGenome=includeGenome, seqChoice=seqChoice, download=download, thorough=thorough)
			if download:
				output.append(sequence)
			else:
				output.append(download is SeqRecord)
		return output

class TerminationPipe(object):
	#Background process class
	def __init__(self, cmd, timeout):
		self.cmd = cmd
		self.timeout = timeout
		self.process = None
		self.output = None
		self.failure = False
	
	def run(self):
		def target():
			self.process = subprocess.Popen(self.cmd, shell=True, stdout=subprocess.PIPE)
			self.output=self.process.communicate()
		
		thread = threading.Thread(target=target)
		thread.start()
		thread.join(self.timeout)
		if thread.is_alive():
			self.process.terminate()
			thread.join()
			self.failure = True

def argsCheck(arguments, parameter, argSplit='-', paramSplit=' '):
	#Check a given string of arguments for a particular parameter, then return its value
	#Assumes my 'standard' "-ARGUMENT VALUE-" setup by kind-of-default
	parameters = arguments.split(argSplit)
	for each in parameters:
		if each == parameter:
			return each.split(parameter+paramSplit)[1]
		else:
			raise RuntimeError("A match value for '" + paramter + "' was not found in the call '" + arguments + "'")

def alignSequences(seqs, method='muscle', tempStem='temp', timeout=99999999):
	#Align given sequences
	#NOTE: Uses subprocess class (above) because internal BioPython methods can hang if you ask for the alignment too soon.
	#NOTE: Requires phylgoeny program to be accessible from everywhere
	#NOTE: beware running this multiple times with the same stem name
	#INPUT: list of SeqRecord objects, method for alignment, temporary stem filename, timeout for pipe (seconds)
	#OUTPUT: alignment object, -or- null tuple if out of time
	#TEST: t = sequenceDownload("vulpes", noSeqs=5, geneName='COI')
	#TEST: alignSequences(t)
	#TEST: alignSequences(t, method='mafft')
	#TEST: alignSequences(t = sequenceDownload("Homo sapiens", noSeqs=10), timeout=5)
	if method == 'muscle':
		inputFile = tempStem + '.fasta'
		outputFile = tempStem + 'Out.fasta'
		commandLine = 'muscle -in ' + inputFile + " -out " + outputFile
		SeqIO.write(seqs, inputFile, "fasta")
	elif method == 'mafft':
		inputFile = tempStem + '.fasta'
		outputFile = tempStem + 'Out.fasta'
		commandLine = 'mafft --auto ' + inputFile + " > " + outputFile
		SeqIO.write(seqs, inputFile, "fasta")
	elif method == 'clustalo':
		inputFile = tempStem + '.fasta'
		outputFile = tempStem + 'Out.fasta'
		commandLine = 'clustalo -i ' + inputFile + " -o " + outputFile + " -v"
		SeqIO.write(seqs, inputFile, "fasta")
	else:
		raise RuntimeError("I haven't implemented anything other than default MUSCLE and mafft alignment methods yet")
	pipe = TerminationPipe(commandLine, timeout)
	pipe.run()
	os.remove(inputFile)
	if not pipe.failure:
		alignment = AlignIO.read(outputFile, 'fasta')
		os.remove(outputFile)
		return alignment
	else:
		raise RuntimeError("Alignment not complete in time allowed")
		return ()

def checkAlignmentList(alignList, method='length', gapType='-'):
	#Check a list of alignments for internal similarity
	def alignLength(alignList):
		return [x.get_alignment_length() for x in alignList]
	
	def noGaps(alignList, gapType):
		meanGapNumber = []
		medianGapNumber = []
		sdGapNumber = []
		maxGapNumber = []
		minGapNumber = []
		for align in alignList:
			gapLength = [len(x) for x in align]
			unGapLength = [len(x.seq.ungap(gapType)) for x in align]
			gapNumber = []
			for g, u in zip(gapLength, unGapLength):
				gapNumber.append(g-u)
			meanGapNumber.append(np.mean(gapNumber))
			medianGapNumber.append(np.median(gapNumber))
			sdGapNumber.append(np.std(gapNumber))
			maxGapNumber.append(max(gapNumber))
			minGapNumber.append(min(gapNumber))
		return {'meanGapNumber':meanGapNumber, 'medianGapNumber':medianGapNumber, 'sdGapNumber':sdGapNumber, 'maxGapNumber':maxGapNumber, 'minGapNumber':minGapNumber}
	
	def gapFraction(alignList, gapType):
		lengths = alignLength(alignList)
		gaps = noGaps(alignList, gapType)
		meanFraction = []
		medianFraction = []
		maxFraction = []
		minFraction = []
		for i in range(len(lengths)):
			meanFraction.append(gaps['meanGapNumber'][i] / lengths[i])
			medianFraction.append(gaps['medianGapNumber'][i] / lengths[i])
			maxFraction.append(gaps['maxGapNumber'][i] / float(lengths[i]))
			minFraction.append(gaps['minGapNumber'][i] / float(lengths[i]))
		return {'meanFraction':meanFraction, 'medianFraction':medianFraction, 'sdGaps':gaps['sdGapNumber'], 'maxFraction':maxFraction, 'minFraction':minFraction}
	
	if method == 'length':
		return alignLength(alignList)
	elif method == 'gapNumber':
		return noGaps(alignList, gapType)
	elif method =='gapFraction':
		return gapFraction(alignList, gapType)
	elif method == 'everything':
		return {'length':alignLength(alignList), 'noGaps':noGaps(alignList, gapType), 'gapFraction':gapFraction(alignList, gapType)}
	else:
		raise RuntimeError("No valid alignment checking method requested")

def checkASequenceList(seqList, method='length', tolerance=None):
	#Check a list of sequences for internal similarity
	def seqLength(seqList):
		return [len(x) for x in seqList]
	
	def quantiles(seqList):
		seqLengths = seqLength(seqList)
		quantileLengths = sp.quantile(seqLengths, [0.05, 0.25, 0.5, 0.75, 0.95])
		maxLength = max(seqLengths)
		minLength = min(seqLengths)
		upperQuantile = [False] * len(seqList)
		lowerQuantile = [False] * len(seqList)
		for i in range(len(seqLengths)):
			if seqLength[i] >= quantileLengths[0]:
				lowerQuantile[i] = True
			elif seqLength[i] <= quantileLengths[4]:
				upperQuantile[i] = True
		return {'seqLengths':lengths, 'maxLength':maxLength, 'minLength':minLength, 'quantiles':quantileLengths, 'lowerQuantile':lowerQuantile, 'upperQuantile':upperQuantile}
	
	if method == 'length':
		return seqLength(seqList)
	elif method == 'quantiles':
		return quantiles(seqList)
	elif method =='quantileDetect':
		if not tolerance:
			raise RuntimeError("Tolerance level not passed to 'checkASequenceList'")
		output = quantiles(seqList)
		upperLimit = output['maxLength'] - output['quantileLengths'][2]
		lowerLimit = output['quantileLengths'][2] - output['minLength']
		if upperLimit > tolerance or lowerLimit > tolerance:
			output['tolerable'] = True
		else:
			output['tolerable'] = False
		return output
	else:
		raise RunetimeError("No valid alignment checking method requested")

def sequenceDisplay(seqList, speciesNames, seqDetect=None):
	#Only works with one gene at the moment
	if(seqDetect):
		maxInputName = max(len(speciesNames))
		print "\nPrinting sequence info. Refer to manual for more details, but note that '^^^' and '___' denote sequences in the upper or lower 5% of sequence lengths'\nGeneral summary:\n"
		if seqDetect['tolerable']:
			print "Sequence lengths within specified tolerance"
		else:
			print "Sequence lengths *outside* specified tolerance"
		print "\nSequence summary:\nInput name".ljust(maxInputName), "GenBankName".ljust(maxSeqName), "Sequence Length".ljust(6), "Issue"
		for i in range(len(seqList)):
			print speciesNames[i].ljust(maxInputName), seqList[i].name.ljust(maxSeqName), len(seqList[i].ljust(6))
	return ()

def phyloGen(alignment, method='RAxML', tempStem='temp', outgroup=None, timeout=None, cladeList=None,  DNAmodel='GTR+G', cleanup=True):
	#Make a phylogeny from a given set of sequences in the background.
	#NOTE: Uses subprocess class (above) because internal BioPython methos can hang if you ask for the alignment too soon.
	#NOTE: Requires phylgoeny program to be accessible from everywhere
	#NOTE: beware running this multiple times with the same stem name, especially with RAxML as it will fail!
	#INPUT: alignment, method for construction, temporary stem filename, timeout for pipe (seconds), whether to clean output files, timeout (if 'None') then everything's prepared for later execution
	#OUTPUT: phylogeny object, -or- null tuple if out of time, -or- a string to type at the command line
	#TEST: t = alignSequences(sequenceDownload("quercus", noSeqs=10, geneName='COI'))
	#TEST: phyloGen(t)
	#TEST: phyloGen(t, method="BEAST")
	#TO-DO: different types of RAxML executable?
	#TO-DO: implement BEAST
	#TO-DO: output some diagnostics from the program while it's running
	################################
	#RAXML##########################
	################################
	if 'RAxML' in method:
		print 'a'
		#RAxML version
		if 'SSE3' in method and 'PTHREADS' in method:
			raxmlVersion = 'raxmlHPC-PTHREADS-SSE3'
		elif 'SSE3' in method:
			raxmlVersion = 'raxmlHPC-SSE3'
		elif 'PTHREADS' in method:
			raxmlVersion = 'raxmlHPC-PTHREADS'
			noThreads = argsCheck(options, 'PTHREADS')
			options += ' -N ' + noThreads
		elif 'localVersion' in method:
			raxmlVersion = 'raxml'
		else:
			raxmlVersion = 'raxmlHPC'
		#DNA model
		if 'noOptimisedModel' in method and 'invariant' in method:
			DNAmodel = ' -m GTRGAMMAI'
		elif 'noOptimisedModel' in method:
			DNAmodel = ' -m GTRGAMMA'
		elif 'invariant' in method:
			DNAmodel = ' -m GTRCATI'
		else:
			DNAmodel = ' -m GTRCAT'
		#Algorithm
		if 'intergratedBootstrap' in method:
			algorithm = ' -f a'
			options += ' -b $RANDOM '
		else:
			algorithm = ' -f d'
		inputFile = tempStem + 'In.phylip'
		outputFile = tempStem + 'Out'
		fileLine = ' -s ' + inputFile + ' -n ' + outputFile
		options = ' -p $RANDOM'
		#Outgroup(s), assuming they're in the right format for RAxML!
		if outgroup:
			options += ' -o ' + outgroup
		AlignIO.write(alignment, inputFile, "phylip")
		commandLine = raxmlVersion + fileLine + algorithm + DNAmodel + options
		if not timeout:
			return commandLine
	elif 'BEAST' in method:
		################################
		#BEAST##########################
		################################
		def indent(elem, level=0):
		    i = "\n" + level*"  "
		    if len(elem):
		        if not elem.text or not elem.text.strip():
		            elem.text = i + "  "
		        if not elem.tail or not elem.tail.strip():
		            elem.tail = i
		        for elem in elem:
		            indent(elem, level+1)
		        if not elem.tail or not elem.tail.strip():
		            elem.tail = i
		    else:
		        if level and (not elem.tail or not elem.tail.strip()):
		            elem.tail = i
		
		def createElement(name, id):
			element=ET.Element("taxa", id=cladeNames[i])
			element.tail="\n"
			return element
		
		if DNAmodel == 'GTR+G':
			#Load the base XML - DIRTY HACK!
			baseXML = ET.parse("/Users/will/Documents/code/phylogen/trunk/GTRBEAST.xml")
			tree = XML.getroot()
			#Add each clade (often one), and its species, with stupid clade names
			for clade in range(len(cladeList)):
				tempClade = createElement('taxa', str(i))
				for sp in clade:
					tempSp = createElement('taxon', sp)
					tempClade.append(tempSp)
				tree.insert(1, temp)
			#Make the bulk of the XML
			cpJumper=True
			coalJumper=True
			priorJumper=True
			logJumper=True
			for i in range(len(tree)):
				if tree[i].tag=="coalescentTree":
					if coalJumper:
						coalJumper=False
						tree[i].remove(tree[i][0])
						temp=ET.Element("constrainedTaxa")
						ET.SubElement(temp, "taxa", idref="taxa")
						for clade in cladeNames:
							subTemp=ET.Element("tmrca", monophyletic="true")
							ET.SubElement(subTemp, "taxa", idref=clade)
							temp.append(subTemp)
							tree[i].append(temp)
				if tree[i].tag=="compoundParameter":
					if cpJumper:
						cpJumper=False
						for clade in cladeNames:
							temp=ET.Element("tmrcaStatistic", id="tmrca_" + clade, includeStem="false")
							subTemp=ET.Element("mrca")
							ET.SubElement(subTemp, "taxa", idref=clade)
							temp.append(subTemp)
							ET.SubElement(temp, "treeModel", idref="treeModel")
							tree.insert(i, temp)
						for clade in cladeNames:
							temp=ET.Element("monophylyStatistic", id="monophyly_" + clade)
							subTemp=Element("mrca")
							ET.SubElement(subTemp, "taxa", idref=clade)
							temp.append(subTemp)
							ET.SubElement(temp, "treeModel", idref="treeModel")
							tree.insert(i, temp)
				if tree[i].tag=="prior":
					if priorJumper:
						priorJumper=False
						temp=ET.Element("booleanLikelihood")
						for clade in cladeNames:
							ET.SubElement(temp, "monophylyStatistic", idref="monophyly_"+clade)
						tree[i].append(temp)
				if tree[i].tag=="logTree":
					if logJumper:
						logJumper=False
						for clade in cladeNames:
							ET.SubElement(tree[i], "tmrcaStatistic", idref="tmarca_"+clade)
			indent(tree)
			ET.ElementTree(tree).write("/Users/will/Documents/rbcL Phylogeny/sequences5Test.xml")
	else:
		raise RuntimeError("Construction method must be RAxML")
	pipe = TerminationPipe(commandLine, timeout)
	pipe.run()
	if not pipe.failure:
		tree = Phylo.read('RAxML_bestTree.' + outputFile, "newick")
		if cleanup:
			os.remove(inputFile)
			dirList = os.listdir(os.getcwd())
			for each in dirList:
				if re.search("(RAxML)", each):
					os.remove(each)
		return tree
	else:
		raise RuntimeError("Either phylogeny building program failed, or ran out of time")

def trimSequence(seq, DNAtype='Standard', gapType='-'):
	#Trim a sequence to contain only a coding region, and make it start on the open reading frame
	#NOTE: Assumes the longest fragment is the ORF (as does GenBank)
	#INPUT: sequence, DNA type
	#OUTPUT: SeqRecord
	#TEST: t = alignSequences(sequenceDownload("quercus", noSeqs=10, geneName='COI'))
	#TEST: trimSequence(sequenceDownload("quercus robur"))
	#TO-DO: this is VERY clunky code; a priority to neaten
	stop = CodonTable.unambiguous_dna_by_name[DNAtype].stop_codons
	start = CodonTable.unambiguous_dna_by_name[DNAtype].start_codons
	directionFrame = []
	output = str()
	nexus = []
	for frame, position in zip(['f', 'r'] * 3, [0, 1, 2] * 2):
		holder = seq.seq.reverse_complement().tostring()[position : len(seq)] if frame=='r' else seq.seq.tostring()[position : len(seq)]
		k = 0
		temp = ['']
		startLock = False
		while len(holder) > 2:
			if startLock:
				temp[k] = temp[k] + holder[0:3]
				if holder[0:3] in stop:
					temp.append('')
					k = k + 1
					startLock = False
			else:
				if holder[0:3] in start:
					temp[k] = temp[k] + holder[0:3]
					startLock = True
			holder = holder[3 : len(holder)]
		temp = sorted(temp, key=len, reverse=True)
		if len(temp[0]) > len(output):
			output = temp[0]
			directionFrame = frame + str(position)
	seq = SeqRecord(Seq(output), id=seq.id, description=seq.description)
	return(seq)

def findGeneInSeq(seq, gene, trimSeq=False, DNAtype='Standard', gapType='-'):
	#Search for a given gene in a given sequence, and return it (or raise an error)
	#REQUIRES the sequence be in GenBank (gb) format, or there aren't any features to search!
	#Can optionally trim the sequence if you wish
	#NOTE: COI can often be recorded as COX1...
	if seq.features:
		for feature in seq.features:
			if 'gene' in feature.qualifiers.keys():
				if gene in feature.qualifiers['gene']:
					extractor = SeqFeature(feature.location)
					foundSeq = extractor.extract(seq)
					if trimSeq:
						return trimSequence(foundSeq, DNAtype=DNAtype, gapType=gapType)
					else:
						return foundSeq
		else:
			raise RuntimeError('Gene not found in sequence')
	else:
		raise RuntimeError('No sequence features found: are you using a GenBank record?')

def rateSmooth(phylo, method='PATHd8', nodes=tuple(), sequenceLength=int(), tempStem='temp', timeout=999999):
	#Rate smooth a given phylogeny
	#NOTE: assuming (checked...ish) that getTipNames always returns the outgroup as the last species name
	#INPUT: phylogeny, method, sequence length, (future) nodes over which to smooth and their ages
	#OUTPUT: smoothed phylogeny -or- error
	#TEST: rateSmooth(tree, sequenceLength=1406) # will get some errors as that's not likely to be the right sequence length...
	#TEST: trimSequence(sequenceDownload("quercus robur"))
	#TO-DO: implement something other than just root-node smoothing (...low priority)
	#TO-DO: neaten up getTipNames
	#TO-DO: quickest way to execute?
	def getTipNames(tree):
		temp = map(lambda x: x.name, tree.get_terminals())
		return([i[0] for i in [each.split("_") for each in temp]])
	
	if not phylo.rooted:
		raise RuntimeError("Phylogeny *must* be rooted")
	if method == 'PATHd8':
		if not nodes:
			if sequenceLength:
				tempPhyloFile = tempStem + 'Phylo'
				phyloText = str()
				Phylo.write(phylo, tempPhyloFile, 'newick')
				with open(tempPhyloFile, 'r') as tempFile:
					phyloText = tempFile.read()
				#phyloText += ';'
				species = getTipNames(phylo)
				mrcaText = '\n\nmrca: ' + species[0] + ', ' + species[len(species)-1] + ', fixage=1;'
				outfile = 'Sequence length = ' + str(sequenceLength) + ';\n\n' + phyloText + mrcaText
				tempPATHd8Input = tempStem + 'PATHd8Input'
				tempPATHd8Output = tempStem + 'PATHd8Output'
				with open(tempPATHd8Input, 'w') as tempFile:
					tempFile.write(outfile)
				commandLine = ' '.join(['PATHd8', tempPATHd8Input, tempPATHd8Output])
				pipe = TerminationPipe(commandLine, timeout)
				pipe.run()
				os.remove(tempPhyloFile)
				os.remove(tempPATHd8Input)
				if not pipe.failure:
					with open(tempPATHd8Output, 'r') as tempFile:
						for i in range(12):
							tempFile.next()
						tempText = tempFile.next()
						tempText = tempText[13:len(tempText)]
						tempSmoothedOutput = tempStem + 'SmoothOutput'
						with open(tempSmoothedOutput, 'w') as tempFile:
							tempFile.write(tempText)
						datedPhylo = Phylo.read(tempSmoothedOutput, 'newick')
					dirList = os.listdir(os.getcwd())
					for each in dirList:
						if re.search('(' + tempStem + ')', each):
							os.remove(each)
					return datedPhylo
			else:
				raise RuntimeError("*Must* provide sequence length")
		else:
			raise RuntimeError("I haven't implemented anything other than defaultPATHd8 smoothing methods yet")
	else:
		raise RuntimeError("I haven't implemented anything other than default PATHd8 smoothing methods yet")


class PhyloGenerator:
	def __init__(self):
		self.fastaFile = ''
		self.GenBankFile = ''
		self.sequences = []
		self.speciesNames = []
		self.downloadInterval = 10
	
	def loadDNAFile(self):
		locker = True
		print "\nIf you have already-downloaded DNA in a single FASTA file, please enter the filename. Otherwise, hit enter to abort"
		while locker:
			inputFile = raw_input("")
			if inputFile:
				try:
					self.sequences.append(list(SeqIO.parse(inputFile, 'fasta')))
					self.fastaFile = inputFile
					locker = False
				except IOError:
					print "\nFile not found. Try again, or hit end-of-file to exit"
			else:
				print "\nNo DNA loaded"
				locker = False
	
	def loadGenBank(self):
		locker = True
		aborted = False
		print "\nIf you need to download gene sequences from GenBank, please enter the filename of the species list (each species on a new line). Otherwise, hit enter to abort"
		while locker:
			inputFile = raw_input("")
			if inputFile:
				try:
					with open(inputFile, 'r') as f:
						for each in f:
							self.speciesNames.append(each)
					
					self.genbankFile = inputFile
					locker = False
				except IOError:
					print "\nFile not found. Try again, or hit end-of-file to exit"
			else:
				print "\nNo DNA downloaded"
				locker = False
				aborted = True
		if not aborted:
			print "\n", len(self.speciesNames), "Species loaded."
			print "\nPlease enter a valid email address to let Entrez know who you are. It's *your* fault if this is not valid, and you will likely have your IP address barred from using GenBank if you don't enter one"
			self.EntrezEmail = raw_input("")
			print"\nPlease enter the name of the gene you want to use, e.g. 'COI' for cytochrome oxidase one, or just hit enter to abort"
			inputGene = raw_input("")
			if inputGene:
				locker = False
				for species in self.speciesNames:
					temp = sequenceDownload(species, inputGene)
					self.sequences.append(temp)
					time.sleep(self.downloadInterval)
			else:
				self.GenBankFile = ''
				self.speciesNames = []
				print "\nGenBank Download aborted"
	
	def DNALoaded():
		if self.sequences:
			pass
		else:
			print "\nNo DNA loaded. Exiting."
			sys.exit()
	
	def dnaChecking():
		pass
	

def main():
	currentState = PhyloGenerator()
	print "\n\nWelcome to phyloGenerator! Let's make a phylogeny!"
	print "\nDNA INPUT"
	currentState.loadDNAFile()
	print "\nDNA DOWNLOAD"
	currentState.loadGenBank()
	currentState.DNALoaded()
	print "\nDNA CHECKING"
	currentState.dnaChecking()

if __name__ == '__main__':
	main()
