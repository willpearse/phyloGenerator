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
import scipy as sp #Stats (quantiles, etc.)
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
		if verbose: progBar = pbar.ProgressBar(widgets=[pbar.Percentage(), pbar.Bar()], maxval=len(speciesList)).start()
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
	
	def run(self, silent=False):
		def target():
			if silent:
				self.process = subprocess.Popen(self.cmd, shell=True)
			else:
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

def alignSequences(seqs, method='muscle', tempStem='temp', timeout=99999999, silent=False):
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
		raise RuntimeError("Alignment method must be 'muscle', 'mafft' or 'clustalo'.")
	pipe = TerminationPipe(commandLine, timeout)
	pipe.run(silent=silent)
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
		return {'mean':meanGapNumber, 'median':medianGapNumber, 'sd':sdGapNumber, 'max':maxGapNumber, 'min':minGapNumber}
	
	def gapFraction(alignList, gapType):
		lengths = [x.get_alignment_length() for x in alignList]
		gaps = noGaps(alignList, gapType)
		meanFraction = []
		medianFraction = []
		maxFraction = []
		minFraction = []
		sdFraction = []
		for i in range(len(lengths)):
			meanFraction.append(gaps['mean'][i] / float(lengths[i]))
			medianFraction.append(gaps['median'][i] / float(lengths[i]))
			maxFraction.append(gaps['max'][i] / float(lengths[i]))
			minFraction.append(gaps['min'][i] / float(lengths[i]))
		return {'mean':meanFraction, 'median':medianFraction, 'max':maxFraction, 'min':minFraction}
	
	if method == 'length':
		return alignList.get_alignment_length()
	elif method == 'gapNumber':
		return noGaps(alignList, gapType)
	elif method =='gapFraction':
		return gapFraction(alignList, gapType)
	elif method == 'everything':
		return {'length':[x.get_alignment_length() for x in alignList], 'noGaps':noGaps(alignList, gapType), 'gapFraction':gapFraction(alignList, gapType)}
	else:
		raise RuntimeError("No valid alignment checking method requested")

def checkSequenceList(seqList, method='length', tolerance=None):
	#Check a list of sequences for internal similarity
	def seqLength(seqList):
		return [len(x) for x in seqList]
	
	def quantiles(seqList):
		seqLengths = seqLength(seqList)
		quantileLengths = sp.percentile(seqLengths, [5, 25, 5, 75, 95])
		maxLength = max(seqLengths)
		minLength = min(seqLengths)
		upperQuantile = [False] * len(seqList)
		lowerQuantile = [False] * len(seqList)
		for i in range(len(seqLengths)):
			if seqLengths[i] <= quantileLengths[0]:
				lowerQuantile[i] = True
			elif seqLengths[i] >= quantileLengths[4]:
				upperQuantile[i] = True
		return {'seqLengths':seqLengths, 'maxLength':maxLength, 'minLength':minLength, 'quantileLengths':quantileLengths, 'lowerQuantile':lowerQuantile, 'upperQuantile':upperQuantile}
	
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
		#Get longest species name
		maxInputName = 0
		for each in speciesNames:
			if each:
				if len(each) > maxInputName:
					maxInputName = len(each)
		#Get longest sequence name
		maxSeqName = 0
		for each in seqList:
			if each:
			 	if len(each.name) > maxSeqName:
					maxSeqName = len(each.name)
		if maxSeqName < len("GenBankName"): maxSeqName = len("GenBankName")
		#Print out details, doing things differently if we haven't analysed them
		if(seqDetect):
			print "\nPrinting sequence info. Refer to manual for more details, but note that '^^^' and '___' denote sequences in the upper or lower 5% of sequence lengths'\nGeneral summary:\n"
			if seqDetect['tolerable']:
				print "Sequence lengths within specified tolerance"
			else:
				print "Sequence lengths *outside* specified tolerance"
			print "\nSequence summary:\nSeq ID ".ljust(len("Seq ID ")), "Input name".ljust(maxInputName), "GenBankName".ljust(maxSeqName), "Sequence Length".ljust(6)
			for i in range(len(seqList)):
				if seqList[i]:
					if seqDetect['upperQuantile'][i]:
						stars = "^^^"
					elif seqDetect['lowerQuantile'][i]:
						stars = "___"
					else:
						stars = ""
					print str(i).ljust(len("Seq ID ")), str(speciesNames[i]).ljust(maxInputName), str(seqList[i].name).ljust(maxSeqName), str(len(seqList[i])).ljust(6), stars
				else:
					print str(i).ljust(len("Seq ID ")), str(speciesNames[i]).ljust(maxInputName), "NO SEQUENCE".ljust(maxSeqName)
		else:
			print "\nPrinting sequence info. Refer to manual for more details."
			print "\nSequence summary:\nSeq ID ".ljust(len("Seq ID ")), "Input name".ljust(maxInputName), "GenBankName".ljust(maxSeqName), "Sequence Length".ljust(6)
			for i in range(len(seqList)):
				if seqList[i]:
					print str(i).ljust(len("Seq ID ")), str(speciesNames[i]).ljust(maxInputName), str(seqList[i].name).ljust(maxSeqName), str(len(seqList[i])).ljust(6)
				else:
					print str(i).ljust(len("Seq ID ")), str(speciesNames[i]).ljust(maxInputName), "NO SEQUENCE".ljust(maxSeqName), "NA".ljust(6)

def alignmentDisplay(alignList, alignMethods, alignDetect=None):
	assert len(alignList)==len(alignMethods)
	print "Below are details of the alignment. Please refer to the manual for more details."
	if alignDetect:
		print "ID", "Alignment".ljust(12), "Length".ljust(10), "Med. Gaps".ljust(15), "SD Gaps".ljust(15), "Min-Max Gaps".ljust(15), "Med. Gap Frac.".ljust(15), "Max Gap Frac.".ljust(15)
		for i in range(len(alignList)):
			print '{ID:3}{alignment:<12}{length:<10}{medGaps:<16.1f}{sdGaps:<16.2f}{minMaxGaps:15}{medGapsFrac:<16.2f}{minMaxGapsFrac:15}'.format(ID=str(i), alignment=alignMethods[i], length=alignDetect['length'][i], medGaps=alignDetect['noGaps']['median'][i], sdGaps=alignDetect['noGaps']['sd'][i], minMaxGaps=str(str(round(alignDetect['noGaps']['min'][i],3))+" - "+str(round(alignDetect['noGaps']['max'][i],3))), medGapsFrac=alignDetect['gapFraction']['median'][i], minMaxGapsFrac=str(str(round(alignDetect['gapFraction']['min'][i],3))+"-"+str(round(alignDetect['gapFraction']['max'][i],3))))
	else:
		print "ID", "Alignment", "Length"
		for i in range(len(alignList)):
			print str(i).ljust(len("ID")), alignMethods[i].ljust(len("Alignment")), str(alignList[i].get_alignment_length()).ljust("Length")

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
		#RAxML version
		if 'SSE3' in method and 'PTHREADS' in method:
			raxmlCompile = '-PTHREADS-SSE3'
		elif 'SSE3' in method:
			raxmlCompile = '-SSE3'
		elif 'PTHREADS' in method:
			raxmlCompile = '-PTHREADS'
			noThreads = argsCheck(options, 'PTHREADS')
			options += ' -N ' + noThreads
		else:
			raxmlCompile = ''
		if 'localVersion' in method:
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
		commandLine = raxmlVersion + raxmlCompile + fileLine + algorithm + DNAmodel + options
		if not timeout:
			return commandLine
	elif 'BEAST' in method:
		################################
		#BEAST##########################
		################################
		def indent(elem, level=0):
			i = "\n" + level*"	"
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
		if 'RAxML' in method:
			tree = Phylo.read('RAxML_bestTree.' + outputFile, "newick")
			if cleanup:
				os.remove(inputFile)
				dirList = os.listdir(os.getcwd())
				for each in dirList:
					if re.search("(RAxML)", each):
						os.remove(each)
					if tempStem+"In.phylip.reduced"==each:
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
	def __init__(self, stem):
		self.fastaFile = ''
		self.GenBankFile = ''
		self.sequences = []
		self.speciesNames = []
		self.downloadInterval = 2
		self.stem = stem
	
	def loadDNAFile(self):
		locker = True
		print "\nIf you have already-downloaded DNA in a single FASTA file, please enter the filename. Otherwise, hit enter to continue."
		while locker:
			inputFile = raw_input("")
			if inputFile:
				try:
					tempSeqs = list(SeqIO.parse(inputFile, 'fasta'))
					self.sequences.extend(tempSeqs)
					self.speciesNames.extend([x.id for x in tempSeqs])
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
							self.speciesNames.append(each.strip())
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
			Entrez.email = raw_input("")
			print"\nPlease enter the name of the gene you want to use, e.g. 'COI' for cytochrome oxidase one, or just hit enter to abort"
			inputGene = raw_input("")
			if inputGene:
				locker = False
				for species in self.speciesNames:
					temp = sequenceDownload(species, inputGene)
					self.sequences.append(temp)
					print "...", species, "found"
					time.sleep(self.downloadInterval)
			else:
				self.GenBankFile = ''
				self.speciesNames = []
				print "\nGenBank Download aborted"
	
	def DNALoaded(self):
		if not self.sequences:
			print "\nNo DNA loaded. Exiting."
			sys.exit()
	
	def dnaChecking(self, tolerance=0.1):
		self.tolerance = tolerance
		self.dnaCheck = checkSequenceList(self.sequences, tolerance=self.tolerance, method="quantileDetect")
		sequenceDisplay(self.sequences, self.speciesNames, self.dnaCheck)
	
	def dnaEditing(self):
		locker = True
		print "\nTo delete a sequence, enter its SeqID and press return - *one sequence at a time*\nTo continue, press enter with no SeqID\n"
		while(locker):
			inputSeq = raw_input("")
			if inputSeq:
				if int(inputSeq) in range(len(self.sequences)):
					del self.sequences[int(inputSeq)]
					print "SeqID", inputSeq, "Successfully deleted"
					print "Re-calulating summary statistics..."
					self.dnaChecking()
				else:
					print "Sorry, I didn't recognise", inputSeq, "- try again"
			else:
				locker = False
				print "No more sequences to delete. Continuing."
	
	def alignmentEditing(self):
		locker = True
		print "\nTo delete a sequence, enter its SeqID and press return - *one sequence at a time*\nTo continue, press enter with no SeqID\n"
		while(locker):
			inputSeq = raw_input("")
			if inputSeq:
				if int(inputSeq) in range(len(self.sequences)):
					del self.sequences[int(inputSeq)]
					print "SeqID", inputSeq, "Successfully deleted"
					print "Re-calulating alignment and summary statistics..."
					self.align()
					self.dnaChecking()
				else:
					print "Sorry, I didn't recognise", inputSeq, "- try again"
			else:
				locker = False
				print "No more sequences to delete. Continuing."
	
	def align(self):
		print "Enter the name of an alignment method ('muscle', 'mafft', 'clustalo'), 'everything' to do all three and compare their outputs, or simply hit return to align with muscle."
		locker = True
		methods = ['muscle', 'mafft', 'clustalo']
		while locker:
			alignInput = raw_input("")
			if alignInput:
				if alignInput in methods:
					print "Aligning DNA with default settings of", alignInput
					self.alignment = alignSequences(self.sequences, method=alignInput, tempStem='temp', timeout=99999999)
					self.alignmentMethod = alignInput
					print "\nAlignment complete!"
					locker = False
				elif alignInput == "everything":
					print "Aligning DNA with:"
					print "...MUSCLE"
					self.alignmentList = [alignSequences(self.sequences, method="muscle", tempStem='temp', timeout=99999999)]
					print "\n...MAFFT"
					self.alignmentList.append(alignSequences(self.sequences, method="mafft", tempStem='temp', timeout=99999999))
					print "\n...Clustal-O"
					self.alignmentList.append(alignSequences(self.sequences, method="clustalo", tempStem='temp', timeout=99999999))
					self.alignmentListNames = ['MUSCLE', 'MAFFT', 'Clustal-O']
					print "\nAlignments complete!"
					self.alignmentChoice()
				else:
					print "Sorry, I didn't recognise", alignInput, "- please try again."
			else:
				print "Alignging DNA with default settings of MUSCLE"
				self.alignment = alignSequences(self.sequences, method="muscle", tempStem='temp', timeout=99999999)
				self.alignmentMethod = "muscle"
				print "\nAlignment complete!"
	
	def alignmentChoice(self):
		self.alignmentCheck = checkAlignmentList(self.alignmentList, method='everything')
		alignmentDisplay(self.alignmentList, self.alignmentListNames, self.alignmentCheck)
		print "\nAny problems with your alignment likely result from poor-quality DNA sequences."
		print "To write out your alignments and view these sequences in an external viewer, type 'output'"
		print "To edit DNA sequences, type 'edit'. To choose an alignment, enter its ID number."
		locker = True
		while locker:
			alignInput = raw_input("Alignment Choice Input: ")
			if alignInput == "output":
				for alignment, method in zip(self.alignmentList, self.alignmentListNames):
					AlignIO.write(alignment, self.stem+"_TEMP_alignment_"+method+".fasta", "fasta")
				print "Alignments written to your working directory."
			elif alignInput == "edit":
				self.dnaEditing()
			elif alignInput in range(self.alignmentList):
				self.alignment = self.alignmentList[alignInput]
				self.alignmentMethod = self.alignmentListNames[alignInput]
				print self.alignmentMethod + "alignment chosen."
				locker = False
	
	def phylogen(self, method="RAxML-localVersion"):
		print"\n Running with options:", method
		self.phylogeny = phyloGen(self.alignment, method=method, timeout=999)
		print"\nRun complete!"
	
	def rateSmooth(self):
		print "\nIf you require a rate-smoothed version of your phylogeny, type the name of the outgroup below.\nOtherwise, hit enter to continue without rate smoothing."
		spNames = [x.name for x in self.phylogeny.get_terminals()]
		locker = True
		while locker:
			inputSmooth = raw_input("")
			if inputSmooth:
				if inputSmooth in spNames:
					self.phylogeny.root_with_outgroup(inputSmooth)
					self.smoothPhylogeny = rateSmooth(self.phylogeny, sequenceLength=self.alignment.get_alignment_length())
					locker = False
			else:
				self.smoothPhylogeny = False
	
	def cleanUpSequences(self):
		cleaned = []
		for i in reversed(range(len(self.sequences))):
			if not self.sequences[i]:
				cleaned.append(self.speciesNames[i])
				del self.sequences[i]
				del self.speciesNames[i]
		if cleaned:
			print "\nThe following species did not have any DNA associated with them, and so have been excluded:"
			for each in cleaned:
				print "\n", each
	
	def renameSequences(self):
		self.genBankIDs = []
		for i in range(len(self.sequences)):
			self.genBankIDs.append(self.sequences[i].id)
			self.sequences[i].name = self.speciesNames
	
	def writeOutput(self):
		#Alignment
		AlignIO.write(self.alignment, self.stem+"_alignment.fasta", 'fasta')
		#Sequence info
		with open(self.stem+"_sequence_info.txt", 'w') as f:
			f.write("Sequence ID, Species Name\n")
			for i in range(len(self.sequences)):
				f.write(self.genBankIDs[i]+","+self.speciesNames[i]+"\n")
		#Phylogeny
		Phylo.write(self.phylogeny, self.stem+"_phylogeny.tre", 'newick')
		if self.smoothPhylogeny:
			Phylo.write(self.smoothPhylogeny, self.stem+"_phylogeny_smoothed.tre", 'newick')
	

def main():
	print "\n\nWelcome to phyloGenerator! Let's make a phylogeny!"
	print "---Please go to http://willpearse.github.com/phyloGenerator for help"
	print "---Written by Will Pearse (will.pearse@gmail.com)"
	print "\nLet's get going!\nPlease input a 'stem' name for all your output (phylogeny, sequences, etc.)"
	stem = raw_input("")
	currentState = PhyloGenerator(stem=stem.strip())
	print "\nDNA INPUT"
	currentState.loadDNAFile()
	print "\nDNA DOWNLOAD"
	currentState.loadGenBank()
	"\nDNA CHECKING"
	currentState.DNALoaded()
	currentState.dnaChecking()
	print "\nYou are now able to delete DNA sequences you have loaded.\nEvery time you delete a sequence, your summary statistics will be re-calculated, and displayed to you again.\n*IMPORTANT*: Sequence IDs may change once you delete a sequence."
	currentState.dnaEditing()
	#TO-DO: allow them to download new sequences for particular species...
	print "\nDNA ALIGNMENT"
	currentState.align()
	print "\nALIGNMENT CHECKING"
	currentState.alignmentChecking()
	print "\nYou are again able to delete DNA sequences you have loaded.\nEvery time you delete a sequence, your alignment and statistics will be re-calculated, and displayed to you again.\n*IMPORTANT*: Sequence IDs may change once you delete a sequence."
	currentState.alignmentEditing()
	currentState.cleanUpSequences()
	currentState.renameSequences()
	print "\nPHYLOGENY GENERATION"
	currentState.phylogen()
	currentState.rateSmooth()
	currentState.writeOutput()
	print "\nCongratulations! Exiting phyloGenerator."

if __name__ == '__main__':
	main()
