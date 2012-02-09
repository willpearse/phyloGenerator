#!/usr/bin/env python
# encoding: utf-8
"""
module.py
First attempt at putting together an automated phylogeny assembly programme for public demonstration.
This version will be funcitonal, i.e. there will be no class definitions
Created by Will Pearse on 2011-08-24.
Copyright (c) 2011 Imperial College london. All rights reserved.
TO-DO:
* Doc-strings
* Unit tests
* Find relatives
* TerminationPipe quiet and better controlled
* BEAST
* 'Execute later' code
* Options file
* Generate constraint tree from GenBank, Phylomatic, web servers, etc.
* Be able to delete a particular gene from the line-up
"""
import pdb
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
import time #For waiting between sequence downloads
import argparse #For command line arguments
import webbrowser #Load website on request
import sys #To exit on errors
import urllib #Download files for install

def taxonIDLookup(taxonID):
	#Lookup a species ID in the NCBI taxonomy database
	#INPUT: species ID
	#OUTPUT: tuple of scientific name, and (tuple) of lineage info (sorted genus-order), taxon ID, mitocondrial code of taxon
	#			-OR- empty tuple
	#NOTE: order of return isn't what you'd want, probably, but done like this for historical reasons
	#NOTE: you wouldn't want to look up anything other than a species with this
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

def findLineage(spName):
	try:
		handleSpName = Entrez.esearch(db="taxonomy", term=spName)
		resultsSpName = Entrez.read(handleSpName)
		handleSpName.close()
		handleID = Entrez.efetch(db="Taxonomy", id=resultsSpName['IdList'], retmode="xml")
		resultsSpID = Entrez.read(handleID)
		lineage = resultsSpID[0]["Lineage"].split("; ")
		lineage.append(spName)
		lineage.reverse()
		return lineage
	except:
		return ()

def findRelativeSequence(spName, geneName=None, cladeDepth=0, thorough=False, rettype='gb', titleText=None, noSeqs=1, seqChoice='random', download=True, retStart=0, retMax=20, targetLength=None, trimSeq=False, DNAtype='Standard', gapType='-', includeGenome=True, includePartial=True):
	#Download a relative's sequence (to be used if one can't be found for that target species)
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
				#Need to check to see if this gene is actually in this sequence (...)
				try:
					seq = findGeneInSeq(seq, geneName, trimSeq=trimSeq, DNAtype=DNAtype, gapType=gapType)
				except:
					seq = dwnSeq(includeGenome=True, includePartial=True)
					if seq:
						return seq
					else:
						return ()
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
	def findBestGene(foundSeqsArray):
		geneHits = foundSeqsArray.sum(axis=0)
		for i in range(len(geneHits)):
			if geneHits[i] == max(geneHits):
				return i
	
	if type(geneNames) is list:
		if targetNoGenes == len(geneNames):
			targetNoGenes = None
		
		#Download number of genes and histories for each species
		searchResults = []
		foundSeqs = []
		foundSeqsBool = np.zeros((len(speciesList), len(geneNames)), int)
		for i in range(len(speciesList)):
			if verbose: print "Searching for:", speciesList[i]
			speciesGenes = []
			for k in range(len(geneNames)):
				sequence = sequenceDownload(speciesList[i], geneNames[k], titleText=titleText, noSeqs=noSeqs, includePartial=includePartial, includeGenome=includeGenome, seqChoice=seqChoice, download=download, thorough=thorough)
				if download:
					speciesGenes.append(sequence)
				foundSeqsBool[i,k] = 1
			if download:
				foundSeqs.append(speciesGenes)
		if targetNoGenes:
			currentFoundSeqs = foundSeqsBool
			currentGeneNames = geneNames[:]
			bestGenes = []
			for each in range(targetNoGenes):
				bestGeneIndex = findBestGene(currentFoundSeqs)
				bestGenes.append(currentGeneNames[bestGeneIndex])
				currentFoundSeqs = np.delete(currentFoundSeqs, bestGeneIndex, axis=1)
				del currentGeneNames[bestGeneIndex]
			output = []
			for i,gene in enumerate(bestGenes):
				currentGene = []
				for j,sp in enumerate(foundSeqs):
					currentGene.append(foundSeqs[i][j])
				output.append(currentGene)
			return (output, bestGenes)
		else:
			return foundSeqsBool
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

def alignSequences(seqList, method='muscle', tempStem='temp', timeout=99999999, silent=False, nGenes=1):
	finalOutput = []
	output = []
	alignedSomething = False
	if method == 'everything': method = 'muscle-mafft-clustalo-prank'
	for i in range(nGenes):
		geneOutput = []
		seqs = [x[i] for x in seqList]
		if 'muscle' in method:
			inputFile = tempStem + '.fasta'
			outputFile = tempStem + 'Out.fasta'
			commandLine = 'muscle -in ' + inputFile + " -out " + outputFile
			SeqIO.write(seqs, inputFile, "fasta")
			pipe = TerminationPipe(commandLine, timeout)
			pipe.run(silent=silent)
			os.remove(inputFile)
			if not pipe.failure:
				geneOutput.append(AlignIO.read(outputFile, 'fasta'))
				os.remove(outputFile)
				alignedSomething = True
			else:
				raise RuntimeError("MUSCLE alignment not complete in time allowed")
		
		if 'mafft' in method:
			inputFile = tempStem + '.fasta'
			outputFile = tempStem + 'Out.fasta'
			commandLine = 'mafft --auto ' + inputFile + " > " + outputFile
			SeqIO.write(seqs, inputFile, "fasta")
			pipe = TerminationPipe(commandLine, timeout)
			pipe.run(silent=silent)
			os.remove(inputFile)
			if not pipe.failure:
				geneOutput.append(AlignIO.read(outputFile, 'fasta'))
				os.remove(outputFile)
				alignedSomething = True
			else:
				raise RuntimeError("Mafft alignment not complete in time allowed")
		
		if 'clustalo' in method:
			inputFile = tempStem + '.fasta'
			outputFile = tempStem + 'Out.fasta'
			commandLine = 'clustalo -i ' + inputFile + " -o " + outputFile + " -v"
			SeqIO.write(seqs, inputFile, "fasta")
			pipe = TerminationPipe(commandLine, timeout)
			pipe.run(silent=silent)
			os.remove(inputFile)
			if not pipe.failure:
				geneOutput.append(AlignIO.read(outputFile, 'fasta'))
				os.remove(outputFile)
				alignedSomething = True
			else:
				raise RuntimeError("Clustal-o alignment not complete in time allowed")
		
		if 'prank' in method:
			inputFile = tempStem + '.fasta'
			outputFile = tempStem + 'Out.fasta'
			commandLine = 'prank -d=' + inputFile + " -o=" + outputFile
			SeqIO.write(seqs, inputFile, "fasta")
			pipe = TerminationPipe(commandLine, timeout)
			pipe.run(silent=silent)
			os.remove(inputFile)
			if not pipe.failure:
				geneOutput.append(AlignIO.read(outputFile+".2.fas", 'fasta'))
				dirList = os.listdir(os.getcwd())
				for each in dirList:
					if re.search("("+outputFile+")", each):
						os.remove(each)
				alignedSomething = True
			else:
				raise RuntimeError("Prank alignment not complete in time allowed")
		
		output.append(geneOutput)
	
	if not alignedSomething:
		raise RuntimeError("Alignment method must be 'muscle', 'mafft', 'clustalo' or 'prank'.")
	return output

def checkAlignmentList(alignList, method='length', gapType='-'):
	def noGaps(alignList, gapType):
		output = {'mean':[], 'median':[], 'sd':[], 'max':[], 'min':[]}
		for align in alignList:
			gapLength = [len(x) for x in align]
			unGapLength = [len(x.seq.ungap(gapType)) for x in align]
			gapNumber = []
			for g, u in zip(gapLength, unGapLength):
				gapNumber.append(g-u)
			output['mean'].append(np.mean(gapNumber))
			output['median'].append(np.median(gapNumber))
			output['sd'].append(np.std(gapNumber))
			output['max'].append(max(gapNumber))
			output['min'].append(min(gapNumber))
		return output
	
	def gapFraction(alignList, gapType):
		output = {'mean':[], 'median':[], 'max':[], 'min':[]}
		lengths = [x.get_alignment_length() for x in alignList]
		gaps = noGaps(alignList, gapType)
		for i in range(len(lengths)):
			output['mean'].append(gaps['mean'][i] / float(lengths[i]))
			output['median'].append(gaps['median'][i] / float(lengths[i]))
			output['max'].append(gaps['max'][i] / float(lengths[i]))
			output['min'].append(gaps['min'][i] / float(lengths[i]))
		return output
	
	def alignLen(geneAlignList):
		output = []
		for method in geneAlignList:
			output.append(method.get_alignment_length())
		return output
	
	if method == 'length':
		return [alignLen(x) for x in alignList]
	elif method == 'gapNumber':
		return [noGaps(x, gapType) for x in alignList]
	elif method =='gapFraction':
		return [gapFraction(x, gapType) for x in alignList]
	elif method == 'everything':
		return {'length':[alignLen(x) for x in alignList], 'noGaps':[noGaps(x, gapType) for x in alignList], 'gapFraction':[gapFraction(x, gapType) for x in alignList]}
	else:
		raise RuntimeError("No valid alignment checking method requested")

def checkSequenceList(seqList, method='length', tolerance=None):
	#Check a list of sequences for internal similarity
	def seqLength(seqList):
		output = []
		for i in range(len(seqList[0])):
			output.append([len(x[i]) for x in seqList])
		return output
	
	def quantiles(seqList):
		seqLengths = seqLength(seqList)
		output = {'seqLengths':[], 'maxLength':[], 'minLength':[], 'quantileLengths':[], 'lowerQuantile':[], 'upperQuantile':[]}
		for i,each in enumerate(seqLengths):
			output['quantileLengths'].append(sp.percentile(each, [5, 25, 5, 75, 95]))
			output['maxLength'].append(max(each))
			output['minLength'].append(min(each))
			output['upperQuantile'].append([False] * len(each))
			output['lowerQuantile'].append([False] * len(each))
			for k in range(len(each)):
				if each[k] <= output['quantileLengths'][0]:
					output['lowerQuantile'][i][k] = True
				elif each[i] >= output['quantileLengths'][4]:
					output['upperQuantile'][i][k] = True
		return output
	
	if method == 'length':
		return seqLength(seqList)
	elif method == 'quantiles':
		return quantiles(seqList)
	elif method =='quantileDetect':
		if not tolerance:
			raise RuntimeError("Tolerance level not passed to 'checkASequenceList'")
		output = quantiles(seqList)
		output['tolerable'] = []
		for i,each in enumerate(seqList):	
			upperLimit = output['maxLength'][i] - output['quantileLengths'][i][2]
			lowerLimit = output['quantileLengths'][i][2] - output['minLength'][i]
			if upperLimit > tolerance or lowerLimit > tolerance:
				output['tolerable'].append(True)
			else:
				output['tolerable'].append(False)
			return output
	else:
		raise RuntimeError("No valid alignment checking method requested")

def sequenceDisplay(seqList, speciesNames, geneNames, seqDetect=None):
	#Get longest species name
	maxInputName = 0
	for each in speciesNames:
		if each:
			if len(each) > maxInputName:
				maxInputName = len(each)
	maxInputName += 1
	#Get longest sequence name
	maxSeqName = 0
	for each in seqList:
		for seq in each:
			if seq:
				if len(seq.name) > maxSeqName:
					maxSeqName = len(seq.name)
	#if maxSeqName < len("GenBankName"): maxSeqName = len("GenBankName")
	maxSeqName += 1
	#Setup geneNames lengths
	geneNameLengths = [len(x) for x in geneNames]
	for i, each in enumerate(geneNameLengths):
		if each < 6:
			geneNameLengths[i] = 6
	
	#Print out details, doing things differently if we haven't analysed them
	if(seqDetect):
		print "\nPrinting sequence info. Refer to manual for more details, but note that '^^^' and '___' denote sequences in the upper or lower 5% of sequence lengths'\nGeneral summary:\n"
		if seqDetect['tolerable']:
			print "Sequence lengths within specified tolerance"
		else:
			print "Sequence lengths *outside* specified tolerance"
		print "\nSequence summary:\n"
		header = "Sp. ID " + "Input name".ljust(maxInputName)
		for i, each in enumerate(geneNames):
			header += each.ljust(geneNameLengths[i]) + "	 "
		print header
		for i in range(len(seqList)):
			stars = []
			for k in range(len(seqList[i])):
				if seqList[i][k]:
					if seqDetect['upperQuantile'][k][i]:
						stars.append("^^^")
					elif seqDetect['lowerQuantile'][k][i]:
						stars.append("___")
					else:
						stars.append("	 ")
				else:
					#What *does* happen to empty sequences?
					stars.append("	 ")
			
			row = str(i).ljust(len("Seq ID ")) + str(speciesNames[i]).ljust(maxInputName)
			for k in range(len(seqList[i])):
				row += str(len(seqList[i][k])).ljust(geneNameLengths[k]) + " " + stars[k] + " "
			print row
	else:
		print "\nPrinting sequence info. Refer to manual for more details."
		print "\nSequence summary:\n"
		header = "Sp. ID " + "Input name".ljust(maxInputName)
		for i, each in enumerate(geneNames):
			header += each.ljust(geneNameLengths[i])
		print header
		for i in range(len(seqList)):
			row = str(i).ljust(len("Seq ID ")) + str(speciesNames[i]).ljust(maxInputName)
			for k in range(len(geneNames)):
				row += str(len(seqList[i][k])).ljust(geneNameLengths[k])
			print row

def alignmentDisplay(alignments, alignMethods, geneNames, alignDetect=None):
	assert len(alignments)==len(geneNames)
	assert len(alignments[0])==len(alignMethods)
	print "Below are details of the alignment(s). Please refer to the manual for more details."
	for alignNo, alignList in enumerate(alignments):
		print "\nGene", geneNames[alignNo]
		if alignDetect:
			print "ID", "Alignment".ljust(12), "Length".ljust(10), "Med. Gaps".ljust(15), "SD Gaps".ljust(15), "Min-Max Gaps".ljust(15), "Med. Gap Frac.".ljust(15), "Max Gap Frac.".ljust(15)
			for i in range(len(alignList)):
				print  '{ID:3}{alignment:<12}{length:<10}{medGaps:<16.1f}{sdGaps:<16.2f}{minMaxGaps:15}{medGapsFrac:<16.2f}{minMaxGapsFrac:15}'.format(ID=str(i), alignment=alignMethods[i], length=alignDetect['length'][alignNo][i], medGaps=alignDetect['noGaps'][alignNo]['median'][i], sdGaps=alignDetect['noGaps'][alignNo]['sd'][i], minMaxGaps=str(str(round(alignDetect['noGaps'][alignNo]['min'][i],3))+" - "+str(round(alignDetect['noGaps'][alignNo]['max'][i],3))), medGapsFrac=alignDetect['gapFraction'][alignNo]['median'][i], minMaxGapsFrac=str(str(round(alignDetect['gapFraction'][alignNo]['min'][i],3))+"-"+str(round(alignDetect['gapFraction'][alignNo]['max'][i],3))))
		else:
			print "ID", "Alignment", "Length"
			for i in range(len(alignList)):
				print str(i).ljust(len("ID")), alignMethods[i].ljust(len("Alignment")), str(alignList[i].get_alignment_length()).ljust(len("Length"))

def phyloGen(alignment, method='RAxML', tempStem='temp', outgroup=None, timeout=None, cladeList=None,  DNAmodel='GTR+G', constraint=None, cleanup=True):
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
		#Constraint
		if constraint:
			if constraint.is_bifurcating():
				options += " -r " + tempStem + "_constraint.tre"
			else:
				options += " -g " + tempStem + "_constraint.tre"
			Phylo.write(constraint, tempStem + "_constraint.tre", 'newick')
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
	seq = SeqRecord(Seq(output), id=seq.id, description=seq.description, name=seq.name)
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

def cleanAlignment(align, method='trimAl-automated', tempStem='temp', timeout=None):
	if 'trimAl' in method:
		options = ""
		if 'automated' in method:
			options += " -automated1"
		fileLine = " -in " + tempStem + "Input.fasta -out " + tempStem + "Output.fasta -fasta"
		trimalVersion = "trimal"
		commandLine = trimalVersion + fileLine + options
		if timeout:
			pipe = TerminationPipe(commandLine, timeout)
			pipe.run()
			if not pipe.failure:
				align = AlignIO.read(tempStem + "Output.fasta", "fasta")
				os.remove(tempStem + "Output.fasta")
				os.remove(tempStem + "Input.fasta")
				return align
			else:
				raise RuntimeError("Either trimAl failed, or it ran out of time")
		else:
			return commandLine
	else:
		raise RuntimeError("Only automated trimAl methods supported at this time.")

def createConstraintTree(spNames, method="phylomaticTaxonomy", fileName='', tempStem='temp'):
	def recursiveTree(lineageList):
		#Make the current depth's elements
		current = [x[1] for x in lineageList]
		uniqueLevels = list(set(current))
		groupedLists = []
		for each in uniqueLevels:
			groupedLists.append([])
		for currentLineage in lineageList:
			for groupedList, uniqueLevel in zip(groupedLists, uniqueLevels):
				if currentLineage[1] == uniqueLevel:
					groupedList.append(currentLineage)
					break
		#Remove the current taxonomic level (don't pass it on in the recursion)
		if len(groupedLists) == 1:
			return "(" + ",".join([x[0] for x in groupedLists[0]]) + ")"
		else:
			for gList in groupedLists:
				for each in gList:
					del each[1]
		
			return "(" + ",".join([recursiveTree(x) for x in groupedLists]) + ")"
	
	if method == "phylomaticTaxonomy":
		phylogenyFile = ' -f ' + fileName
		with open(tempStem + "taxa.txt", 'w') as taxaFileOutput:
			for sp in spNames:
				taxaFileOutput.write(sp + "\n")
		taxaFile = ' -t ' + tempStem + "taxa.txt"
		commandLine = 'phylomatic' + phylogenyFile + taxaFile
		pipe = TerminationPipe(commandLine, timeout)
		pipe.run(silent=silent)
		os.remove(tempStem + "taxa.txt")
		if not pipe.failure:
			geneOutput.append(AlignIO.read(outputFile, 'fasta'))
		else:
			raise RuntimeError("Phylomatic did not run correctly")
		
	elif method == "GenBank":
		print "THIS DOESN'T WORK!!!!!!!!"
		lineages = [findLineage(x) for x in spNames]
		return recursiveTree(lineages)
	else:
		raise RuntimeError("Unrecognised constraint tree creation method specified")

def createConstraintTreeCaps(spNames):
	def recursiveTree(lineageList):
		#Make the current depth's elements
		current = [x[1] for x in lineageList]
		uniqueLevels = list(set(current))
		groupedLists = []
		for each in uniqueLevels:
			groupedLists.append([])
		for currentLineage in lineageList:
			for groupedList, uniqueLevel in zip(groupedLists, uniqueLevels):
				if currentLineage[1] == uniqueLevel:
					groupedList.append(currentLineage)
					break
		#Remove the current taxonomic level (don't pass it on in the recursion)
		if len(groupedLists) == 1:
			return "(" + ",".join([x[0] for x in groupedLists[0]]) + ")"
		else:
			for gList in groupedLists:
				for each in gList:
					del each[1]
		
		return "(" + ",".join([recursiveTree(x) for x in groupedLists]) + ")"
	
	lineages = [findLineage(x) for x in spNames]
	for each in lineages:
		for i in reversed(range(len(each))):
			if not each[i].istitle():
				del each[i]
	return recursiveTree(lineages)

class TerminationPipe(object):
	#Background process class
	def __init__(self, cmd, timeout):
		self.cmd = cmd
		self.timeout = timeout
		self.process = None
		self.output = None
		self.failure = False
		self.stdout = 'EMPTY'
		self.stderr = 'EMPTY'
	
	def run(self, silent=True):
		def silentTarget():
			tStdout  = open('termPipeStdOut.txt', 'w')
			tStderr  = open('termPipeStdErr.txt', 'w')
			self.process = subprocess.Popen(self.cmd, shell=True, stdout=tStdout, stderr=tStderr)
			self.stdout = open('termPipeStdOut.txt', 'r').readlines()
			self.stderr = open('termPipeStdErr.txt', 'r').readlines()
			os.remove('termPipeStdOut.txt')
			os.remove('termPipeStdErr.txt')
			self.output=self.process.communicate()
		
		def loudTarget():
			self.process = subprocess.Popen(self.cmd, shell=True, stdout=subprocess.PIPE)
			self.output=self.process.communicate()
		
		if silent:
			thread = threading.Thread(target=silentTarget)
		else:
			thread = threading.Thread(target=loudTarget)
		thread.start()
		thread.join(self.timeout)
		if thread.is_alive():
			self.process.terminate()
			thread.join()
			self.failure = True

class PhyloGenerator:
	def __init__(self, stem):
		self.fastaFile = False
		self.GenBankFile = False
		self.sequences = []
		self.speciesNames = []
		self.dnaCheck = []
		self.downloadInterval = 2
		self.stem = stem
		self.phylogeny = False
		self.alignment = []
		self.smoothPhylogeny = False
		self.root = False
		self.genBankIDs = []
		self.constraint = False
		self.genes = []
		self.nGenes = -999
		self.maxGenBankDownload = 50
		self.alignmentMethod = []
		self.alignmentMethodChosen = []
		self.email = ''
		self.codonModels = []
	
	def loadDNAFile(self, inputFile=""):
		if inputFile:
			try:
				tempSeqs = list(SeqIO.parse(inputFile, 'fasta'))
				for each in tempSeqs:
					self.sequences.append([each])
				self.speciesNames.extend([x.name for x in tempSeqs])
				self.fastaFile = inputFile
				if not self.genes:
					print "DNA loaded; please enter the name of the gene you're using below"
					self.genes.append(raw_input("Gene name: "))
					self.nGenes = 1
				self.codonModels.append('Standard')
				print "DNA loaded"
			except IOError:
				print "\nDNA sequence file not found. Exiting..."
				sys.exit()
		else:
			locker = True
			print "\nIf you have already-downloaded DNA in a single FASTA file, please enter the filename. Otherwise, hit enter to continue."
			while locker:
				inputFile = raw_input("")
				if inputFile:
					try:
						tempSeqs = list(SeqIO.parse(inputFile, 'fasta'))
						for each in tempSeqs:
							self.sequences.append([each])
						self.speciesNames.extend([x.name for x in tempSeqs])
						self.fastaFile = inputFile
						if not self.genes:
							print "DNA loaded; please enter the name of the gene you're using below"
							self.genes.append(raw_input("Gene name: "))
							self.nGenes = 1
						self.codonModels.append('Standard')
						print "DNA loaded"
						locker = False
					except IOError:
						print "\nFile not found. Please try again!"
				else:
					print "\nNo DNA loaded"
					locker = False
	
	def loadGenBank(self, inputFile=""):
		if inputFile:
			try:
				with open(inputFile, 'r') as f:
					for each in f:
						self.speciesNames.append(each.strip())
						aborted = False
			except IOError:
				print "\nERROR: File not found. Exiting..."
		else:
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
			if not self.email:
				print "\nPlease enter a valid email address to let Entrez know who you are. It's *your* fault if this is not valid, and you will likely have your IP address barred from using GenBank if you don't enter one"
				Entrez.email = raw_input("")
			if not self.genes:
				print"\nPlease enter the name of the gene you want to use, e.g. 'COI' for cytochrome oxidase one. To enter multiple genes, enter each on a separateline. Just hit enter to abort."
				locker = True
				self.nGenes = 0
				while locker:
					inputGene = raw_input("")
					if inputGene:
						self.genes.append(inputGene)
						self.codonModels.append('Standard')
						self.nGenes += 1
					else:
						locker = False
			if self.genes:
				geneOutput = findGenes(self.speciesNames, self.genes, download=True, seqChoice="medianLength", verbose=True, thorough=True, targetNoGenes=self.nGenes)
				self.sequences = geneOutput[0]
				self.genes = geneOutput[1]
	
	def DNALoaded(self):
		if not self.sequences:
			print "\nNo DNA loaded. Exiting."
			sys.exit()
	
	def dnaChecking(self, tolerance=0.1):
		self.tolerance = tolerance
		self.dnaCheck = checkSequenceList(self.sequences, tolerance=self.tolerance, method="quantileDetect")
		sequenceDisplay(self.sequences, self.speciesNames, self.genes, self.dnaCheck)
	
	def dnaEditing(self):
		def deleteMode(firstTime=True):
			if firstTime:
				print "\nYou're in deletion mode. To delete a species, enter its SeqID and press return.\t*One species at a time please!*"
				print "To change to the 'reload', 'trim' or 'replace' modes, simply type their names then hit enter."
			inputSeq = raw_input("DNA Editing (delete): ")
			if inputSeq:
				try:
					if int(inputSeq) in range(len(self.sequences)):
						del self.sequences[int(inputSeq)]
						del self.speciesNames[int(inputSeq)]
						print "SeqID", inputSeq, "Successfully deleted"
						print "Re-calulating summary statistics..."
						self.dnaChecking()
						return 'delete', False
				except:
					pass
				if inputSeq == "trim":
					return "trim", True
				elif inputSeq == "reload":
					return "reload", True
				elif inputSeq == "replace":
					return "replace", True
				else:
					print "Sorry,", inputSeq, "was not recognised. Please try again."
					return 'delete', False
			else:
				return "EXIT", True
		
		def reloadMode(firstTime=True):
			if not self.email:
				print "\nPlease enter a valid email address to let Entrez know who you are. It's *your* fault if this is not valid, and you will likely have your IP address barred from using GenBank if you don't enter one"
				Entrez.email = raw_input("")
			if firstTime:
				print "\nYou're in reload mode. To reload all sequences in a species, enter its SeqID and press return. To reload just one sequence, enter its SeqID and the gene name."
				print "To reload all sequences above or below 900 bp in length (for example), type '>900' or '<900' respectively."
				print "\t(to reload everything thoroughly, type 'EVERYTHING' (note the caps))"
				print "\t*One species/sequence at a time please!*"
				print "To change to the 'delete', 'trim' or 'reload' modes, simply type their names then hit enter."
			inputSeq = raw_input("DNA Editing (reload):")
			if inputSeq:
				try:
					if int(inputSeq) in range(len(self.sequences)):
						for i,each in enumerate(self.sequences[int(inputSeq)]):
							self.sequences[int(inputSeq)][i] = sequenceDownload(self.speciesNames[int(inputSeq)], self.genes[i], thorough=True, retMax=self.maxGenBankDownload, seqChoice='targetLength', targetLength=self.dnaCheck['quantileLengths'][i][2])
						print "Re-calulating summary statistics..."
						self.dnaChecking()
						return 'reload', False
				except:
					pass
				for i,gene in enumerate(self.genes):
					if gene in inputSeq:
						seqID = re.search("[0-9]*", inputSeq).group()
						if seqID and seqID < len(self.sequences):
							print "Reloading SeqID", seqID, "gene", gene
							self.sequences[seqID][i] = sequenceDownload(self.speciesNames[seqID], self.genes[i], thorough=True, retMax=self.maxGenBankDownload, seqChoice='targetLength', targetLength=self.dnaCheck['quantileLengths'][i][2])
							print "Re-calulating summary statistics..."
							self.dnaChecking()
							return 'reload', False
				if ">" in inputSeq:
					threshold = int(inputSeq.replace('>', ''))
					for i,sp in enumerate(self.sequences):
						for j,gene in enumerate(sp):
							if len(self.sequences[i][j]) > threshold:
								self.sequences[i][j] = sequenceDownload(self.speciesNames[i], self.genes[j], thorough=True, retMax=self.maxGenBankDownload, seqChoice='targetLength', targetLength=self.dnaCheck['quantileLengths'][j][2])
					print "Re-calulating summary statistics..."
					self.dnaChecking()
					return 'reload', False
				if "<" in inputSeq:
					threshold = int(inputSeq.replace('<', ''))
					for i,sp in enumerate(self.sequences):
						for j,gene in enumerate(sp):
							if len(self.sequences[i][j]) < threshold:
								self.sequences[i][j] = sequenceDownload(self.speciesNames[j], self.genes[i], thorough=True, retMax=self.maxGenBankDownload, seqChoice='targetLength', targetLength=self.dnaCheck['quantileLengths'][j][2])
					print "Re-calulating summary statistics..."
					self.dnaChecking()
					return 'reload', False
				if inputSeq == "EVERYTHING":
					geneOutput = findGenes(self.speciesNames, self.genes, download=True, seqChoice="medianLength", verbose=True, thorough=True, retMax=self.maxGenBankDownload, targetNoGenes=self.nGenes)
					self.sequences = geneOutput[0]
					self.genes = geneOutput[1]
					print "Re-calulating summary statistics..."
					self.dnaChecking()
					return 'reload', False
				elif inputSeq == "delete":
					return "delete", True
				elif inputSeq == "trim":
					return "trim", True
				elif inputSeq == "replace":
					return "replace", True
				else:
					print "Sorry,", inputSeq, "was not recognised. Please try again."
					return 'reload', False
			else:
				return "EXIT", True
		
		def trimMode(firstTime=True):
			if firstTime:
				print "\nYou're in trim mode. To trim all sequences in a species, enter its SeqID and press return. To trim just one sequence, enter its SeqID and the gene name."
				print "The type of gene you've downloaded may affect the kind of trimming that takes places. To change/review the type of gene, type 'type' and press enter."
				print "\tTo trim all sequences longer than 1000 bp, type '>1000'"
				print"\t(to trim everything, type 'EVERYTHING' (note the caps))"
				print "\t*One species/sequence at a time please!*"
				print "To change to the 'delete', 'reload' or 'replace' modes, simply type their names then hit enter."
			inputSeq = raw_input("DNA Editing (trim):")
			if inputSeq:
				try:
					if int(inputSeq) in range(len(self.sequences)):
						for i,each in enumerate(self.sequences[int(inputSeq)]):
							if self.sequences[int(inputSeq)]:
								self.sequences[int(inputSeq)][i] = trimSequence(self.sequences[int(inputSeq)][i])
						print "Re-calulating summary statistics..."
						self.dnaChecking()
						return 'trim', False
				except:
					pass
				for i,gene in enumerate(self.genes):
					if gene in inputSeq:
						seqID = re.search("[0-9]*", inputSeq).group()
						if seqID and int(seqID) < len(self.sequences):
							print "Reloading SeqID", seqID, "gene", gene
							self.sequences[seqID][i] = trimSequence(self.sequences[inputSeq][i])
							print "Re-calulating summary statistics..."
							self.dnaChecking()
							return 'trim', False
				if ">" in inputSeq:
					threshold = int(inputSeq.replace('>', ''))
					for i,sp in enumerate(self.sequences):
						for j,gene in enumerate(sp):
							if len(self.sequences[i][j]) > threshold:
								self.sequences[i][j] = trimSequence(self.sequences[i][j])
					print "Re-calulating summary statistics..."
					self.dnaChecking()
					return 'trim', False
				if inputSeq == 'type':
					print "Below are the names and IDs of gene-types."
					print "ID", "Gene Type"
					for i,x in enumerate(CodonTable.unambiguous_dna_by_name.keys()):
						print str(i).ljust(2), x
					print "\nIn the prompt below is the name of a gene. Type the ID number of the codon model for that gene you'd like to use. Note that the standard (default) model is usually number 10."
					for i,gene in enumerate(self.genes):
						locker = True
						while locker:
							choice = raw_input(gene+": ")
							if int(choice) in range(len(CodonTable.unambiguous_dna_by_name.keys())):
								self.codonModels[i] = CodonTable.unambiguous_dna_by_name.keys()[int(choice)]
								locker = False
							else:
								print "Sorry, didn't get that. Try again."
					print "DNA types changed!"
					return 'trim', False
				if inputSeq == "delete":
					return "delete", True
				elif inputSeq == "reload":
					return "reload", True
				elif inputSeq == "replace":
					return "replace", True
				elif inputSeq == "EVERYTHING":
					for i,sp in enumerate(self.sequences):
						for j,gene in enumerate(sp):
							if self.sequences[i][j]:
								self.sequences[i][j] = trimSequence(self.sequences[i][j])
					print "Re-calulating summary statistics..."
					self.dnaChecking()
					return 'trim', False
				else:
					print "Sorry,", inputSeq, "was not recognised. Please try again."
					return 'trim', False
			else:
				return "EXIT", True
		
		def replaceMode(firstTime=True):
			if not self.email:
				print "\nPlease enter a valid email address to let Entrez know who you are. It's *your* fault if this is not valid, and you will likely have your IP address barred from using GenBank if you don't enter one"
				Entrez.email = raw_input("")
			if firstTime:
				print "\nYou're in replace mode. To replace a particular species with a congener, simply type its SeqID and press enter."
				print "\tTo replace all species without any sequences with a congener, type 'EVERYTHING' and press enter."
				print "To change to the 'delete', 'trim' or 'reload' modes, simply type their names then hit enter."
			inputSeq = raw_input("DNA Editing (replace):")
			if inputSeq:
				try:
					if int(inputSeq) in range(len(self.sequences)):
						i = int(inputSeq)
						lineage = findLineage(self.speciesNames[i])
						if lineage:
							replacements = cladeSpecies(lineage[1])
							print "...looking for alternatives for", self.speciesNames[i], "in clade", lineage[1]
							for candidate in replacements:
								temp = []
								locker = False
								for j,gene in enumerate(self.genes):
									temp.append(sequenceDownload(candidate[0], gene))
									if temp[j]:
										locker = True
								if locker:
									self.sequences[i] = temp
									print "......alternative found:", candidate[0], "re-caluclating summary statistics..."
									self.speciesNames[i] = candidate[0] +"(" + self.speciesNames[i] + ")"
									self.dnaChecking()
									return 'replace', False
							else:
								print "No alternative found."
								return 'replace', False
						else:
							print "...Cannot find any entry in GenBank with that name."
							if " " in self.speciesNames[i]:
								newName = self.speciesNames[i].split(' ')[0]
								print "......Trying to find an entry for...", newName
								lineage = findLineage(newName[i])
								if lineage:
									replacements = cladeSpecies(lineage[1])
									print "...looking for alternatives for", newName[i], "in clade", lineage[1]
									for candidate in replacements:
										name = 'ERROR'
										locker = False
										for gene in self.genes:
											temp.append(sequenceDownload(candidate[0], gene))
											if temp:
												locker = True
										if locker:
											self.sequences[i] = temp
											self.speciesNames[i] = candidate[0] +"(" + self.speciesNames[i] + ")"
											print "......alternative found:", candidate[0], "re-caluclating summary statistics..."
											self.dnaChecking()
											return 'replace', False
									else:
										print "No alternative found."
										return 'replace', False
								else:
									print "......Cannot find any entry in GenBank with that name."
							else:
								print "......Can't auto-detect a suitable generic name to search for."
								return 'replace', False
							return 'replace', False
				except:
					pass
				if inputSeq == "EVERYTHING":
					for i,sp in enumerate(self.sequences):
						tracker = 0
						for j,gene in enumerate(sp):
							if gene:
								tracker += 1
						if tracker == 0:
							lineage = findLineage(self.speciesNames[i])
							replacements = cladeSpecies(lineage[1])
							print "...looking for alternatives for", self.speciesNames[i], "in clade", lineage[1]
							for candidate in replacements:
								temp = []
								locker = False
								for gene in self.genes:
									temp.append(sequenceDownload(candidate[0], gene))
									if temp:
										locker = True
								if locker:
									self.sequences[i] = temp
									print "......alternative found:", candidate[0]
								break
					print "Re-calulating summary statistics..."
					self.dnaChecking()
					return 'replace', False
				if inputSeq == "delete":
					return "delete", True
				elif inputSeq == "reload":
					return "reload", True
				elif inputSeq == "trim":
					return "trim", True
				else:
					print "Sorry,", inputSeq, "was not recognised. Please try again."
					return 'replace', False
			else:
				return "EXIT", True
		
		locker = True
		firstTime = True
		mode = "delete"
		while(locker):
			if mode == "delete":
				mode, firstTime = deleteMode(firstTime)
			elif mode == "reload":
				mode, firstTime = reloadMode(firstTime)
			elif mode == "trim":
				mode, firstTime = trimMode(firstTime)
			elif mode == "replace":
				mode, firstTime = replaceMode(firstTime)
			else:
				raise RuntimeError("Unrecognised DNA Editing mode!")
			if mode:
				if mode == "EXIT":
					locker = False
	
	def alignmentEditing(self):
		alignmentDisplay(self.alignment, self.alignmentMethod, self.genes, checkAlignmentList(self.alignment, method='everything'))
		print "\nIt's *strongly* recommended that you take a look at your alignment before continuing; the summary statistics above are unlikely to be sufficient to spot big problems!"
		print "To print out your alignments, type 'output'. To return to the DNA editting stage, type 'DNA', and to align the sequences differently type 'align'."
		print "\t, To automatically trim your sequences using trimAl, type 'trim'"
		print "\nYou *cannot* continue without a chosen alignment for each gene you are using. To choose an alignment, just hit enter"
		locker = True
		while locker:
			inputAlign = raw_input("Alignment Checking:")
			if inputAlign:
				if inputAlign == 'output':
					for i,gene in enumerate(self.genes):
						for j,method in enumerate(self.methods):
							AlignIO.write(self.alignments[i][k], self.stem+"_"+gene+"_"+method, 'fasta')
				elif inputAlign == 'DNA':
					self.dnaChecking()
				elif inputAlign == 'align':
					self.align()
				elif inputAlign == 'trim':
					self.alignment = trimAlignment(self.alignment)
					alignmentDisplay(self.alignment)
				else:
					print "Sorry, I don't understand", inputAlign, "- please try again."
			else:
				locker = False
		print "\nIn the prompt below is the name of a gene. Type the ID number of the alignment you'd like to use for that gene below."
		for i,gene in enumerate(self.genes):
			locker = True
			while locker:
				choice = raw_input(gene+": ")
				if int(choice) in range(len(self.genes)):
					self.alignment[i] = self.alignment[i][int(choice)]
					self.alignmentMethodChosen.append(self.alignmentMethod[int(choice)])
					locker = False
				else:
					print "Sorry, didn't get that. Try again."
	
	def align(self):
		print "Enter the name of an alignment method ('muscle', 'mafft', 'clustalo', 'prank'), 'everything' to do all four and compare their outputs, 'quick' to do only the first three (prank can take a while!), or simply hit return to align with muscle."
		locker = True
		methods = ['muscle', 'mafft', 'clustalo', 'prank']
		while locker:
			alignInput = raw_input("")
			if alignInput:
				if alignInput in methods:
					print "Aligning DNA with default settings of", alignInput
					self.alignment = alignSequences(self.sequences, method=alignInput, tempStem='temp', timeout=99999999, nGenes=len(self.genes))
					self.alignmentMethod.append(alignInput)
					print "\nAlignment complete!"
					locker = False
				elif alignInput == "everything":
					print "Aligning DNA with:"
					print "...MUSCLE"
					self.alignment.append(alignSequences(self.sequences, method="muscle", tempStem='temp', timeout=99999999, nGenes=len(self.genes)))
					print "\n...MAFFT"
					self.alignment.append(alignSequences(self.sequences, method="mafft", tempStem='temp', timeout=99999999, nGenes=len(self.genes)))
					print "\n...Clustal-O"
					self.alignment.append(alignSequences(self.sequences, method="clustalo", tempStem='temp', timeout=99999999, nGenes=len(self.genes)))
					print "\n...Prank"
					self.alignment.append(alignSequences(self.sequences, method="prank", tempStem='temp', timeout=99999999, nGenes=len(self.genes)))
					self.alignmentMethods = ['MUSCLE', 'MAFFT', 'Clustal-O', 'Prank']
					print "\nAlignments complete!"
					locker = False
				else:
					print "Sorry, I didn't recognise", alignInput, "- please try again."
			else:
				print "Alignging DNA with default settings of MUSCLE"
				self.alignment = alignSequences(self.sequences, method="muscle", tempStem='temp', timeout=99999999, nGenes=len(self.genes))
				self.alignmentMethod.append("muscle")
				print "\nAlignment complete!"
				locker = False
	
	def phylogen(self, method="RAxML-localVersion"):
		print"\n Running with options:", method
		self.phylogeny = phyloGen(self.alignment, method=method, constraint=self.constraint, timeout=999)
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
					print "Sorry, I couldn't find", inputSmooth, " in your phylogeny - try again"
			else:
				locker = False
	
	def cleanUpSequences(self):
		cleaned = []
		for i,sp in reversed(list(enumerate(self.sequences))):
			foundSequence = False
			for j,seq in reversed(list(enumerate(sp))):
				if seq:
					foundSequence = True
				else:
					self.sequences[i][j] = "-"
			else:
				if not foundSequence:
					cleaned.append(self.speciesNames[i])
					del self.sequences[i]
					del self.speciesNames[i]
		
		if cleaned:
			print "\nThe following species did not have any DNA associated with them, and so have been excluded:"
			for each in cleaned:
				print "\n", each
	
	def renameSequences(self):
		for i in range(len(self.sequences)):
			tGenBankIDs = []
			for k in range(len(self.sequences[i])):
				tGenBankIDs.append(self.sequences[i][k].id)
				self.sequences[i][k].id = self.speciesNames[i].replace(" ", "_")
			self.genBankIDs.append(tGenBankIDs)
	
	def writeOutput(self):
		#This is unlikely to work if called at any point during execution (alignment's structure changes, etc.)
		# - but this reflects a more fundamental problem in the way in which you've structured everything...
		#Log - TO-DO
		#Sequences
		if self.sequences:
			for i,gene in enumerate(self.genes):
				currentGene = []
				for seq in self.sequences:
					currentGene.append(seq[i])
				SeqIO.write(currentGene, self.stem+"_"+"gene"+"raw_sequences.fasta", 'fasta')
		
		#Alignment
		if self.alignment and self.genes:
			for gene,align in zip(self.genes, self.alignment):
				AlignIO.write(align, self.stem+"_"+gene+"_alignment.fasta", 'fasta')
		
		#Sequence info
		if self.genBankIDs:
			for i,gene in enumerate(self.genes):
				with open(self.stem+"_"+self.genes[i]+"_sequence_info.txt", 'w') as f:
					f.write("Species Name, Sequence ID\n")
					for j,name in enumerate(self.speciesNames):
						f.write(name + "_" + self.genBankIDs[j][i] + "\n")
		#Phylogeny
		if self.phylogeny:
			Phylo.write(self.phylogeny, self.stem+"_phylogeny.tre", 'newick')
		#Constraint tree
		if self.constraint:
			Phylo.write(self.constraint, self.stem+"_constraint.tre", 'newick')
		#Smoothed phylogeny
		if self.smoothPhylogeny:
			Phylo.write(self.smoothPhylogeny, self.stem+"_phylogeny_smoothed.tre", 'newick')
	
	def getConstraint(self, fileName=""):
		if not fileName:
			print "It is *stronlgy* advised that you use a constraint tree with this program."
			print "Please input the filename of your constraint tree (in newick format), or press enter to continue without one."
			locker = True
			while locker:
				inputConstraint = raw_input("Constraint Tree: ")
				if inputConstraint:
					try:
						self.constraint = Phylo.read(inputConstraint, 'newick')
						locker = False
					except IOError:
						print "\nFile not found. Please try again!"
					if self.checkConstraint():
						print "Constraint tree loaded!"
						locker = False
					else:
						print "Constraint tree does *not* match the species names you've inputted. Please load another file."
				else:
					print "...No constraint tree loaded"
					locker = False
		else:
			try:
				self.constraint = Phylo.read(fileName, 'newick')
			except IOError:
				print "\nFile not found. Exiting..."
				sys.exit()
			if self.checkConstraint():
				print "Constraint tree loaded!"
			else:
				print "Constraint tree does *not* match the species names you've inputted. Exiting..."
				sys.exit()
	
	def checkConstraint(self):
		if self.constraint:
			tipLabels = [x.name for x in self.constraint.get_terminals()]
			count = 0
			for each in tipLabels:
				if each in self.speciesNames:
					count += 1
			if count == len(tipLabels):
				return True
			else:
				return False
		else:
			return False
	
	def concatenateSequences(self):
		if len(self.genes) > 1:
			tempAlignment = self.alignment[0]
			for i in range(1, len(self.genes)):
				tempAlignment += self.alignment[i]
			self.alignment = tempAlignment
	
	

def main():
	args = parser.parse_args()
	if args.version:
		print "v0.1"
	elif args.manual:
		 webbrowser.open("http://willpearse.github.com/phyloGenerator")
	else:
		print "\n\nWelcome to phyloGenerator! Let's make a phylogeny!"
		print "---Please go to http://willpearse.github.com/phyloGenerator for help"
		print "---Written by Will Pearse (will.pearse@gmail.com)"
		
		#Stem name
		if args.name:
			print "\nLet's get going!\nUsing stem name", args.name, "..."
			currentState = PhyloGenerator(stem=args.name)
		else:
			print "\nLet's get going!\nPlease input a 'stem' name for all your output (phylogeny, sequences, etc.)"
			stem = raw_input("Stem name: ")
			currentState = PhyloGenerator(stem=stem)
		
		#Email
		if args.email:
			Entrez.email = args.email
			currentState.email = True
		
		#Gene
		if args.gene:
			currentState.genes = args.gene.split(',')
		
		#Gene Number
		if args.nGenes:
			currentState.nGenes = int(args.nGenes)
		
		#Handle sequence input
		if args.alignment and args.dna:
			print "\nERROR: Can't handle an alignment and DNA - suggest you manually strip the alignment and merge the files.\nExiting."
			sys.exit()
		#Alignment
		if args.alignment :
			try:
				currentState.alignment = AlignIO.read(args.dna, 'fasta')
				print "Alignment successfully loaded from file", args.alignment, "..."
			except:
				print "\n!!!Cannot load alignment file! Exiting..."
				sys.exit()
		elif args.dna:
			#Raw DNA
			currentState.loadDNAFile(args.dna)
		elif args.species:
			#Species list
			currentState.loadGenBank(args.species)
		else:
			print "\nDNA INPUT"
			currentState.loadDNAFile()
			print "\nDNA DOWNLOAD"
			currentState.loadGenBank()
		"\nDNA CHECKING"
		currentState.DNALoaded()
		for i in currentState.sequences:
			print "!!!!!!!!!!!!!!!!!!!!!!!!"
			for k in i:
				print k
		currentState.dnaChecking()
		print "\nYou are now able to delete DNA sequences you have loaded.\nEvery time you delete a sequence, your summary statistics will be re-calculated, and displayed to you again.\n*IMPORTANT*: Sequence IDs may change once you delete a sequence."
		currentState.dnaEditing()
		currentState.cleanUpSequences()
		currentState.renameSequences()
		#TO-DO: allow them to download new sequences for particular species...
		
		#Constraint tree
		if args.constraint:
			currentState.getConstraint(args.constraint)
		else:
			print "\nCONSTRAINT TREE"
			currentState.getConstraint()
		
		#Alignment
		
		if not args.alignment:
			print "\nDNA ALIGNMENT"
			currentState.align()
			print "\nALIGNMENT CHECKING"
			currentState.alignmentEditing()
			currentState.concatenateSequences()
			print "\nPHYLOGENY GENERATION"
		currentState.phylogen()
		currentState.rateSmooth()
	#except:
	#	print "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	#	print "Sorry, an unhandled exception has occured."
	#	print "I will try to write out whatever has been done so far."
	#	print "If this isn't obviously something you've done, please contact me - will.pearse@gmail.com"
	#	print "Sorry again!"
	currentState.writeOutput()
	print "\nCongratulations! Exiting phyloGenerator."


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="phyloGenerator - phylogeny generation for ecologists.", epilog="Help at http://willpearse.github.com/phyloGenerator - written by Will Pearse")
	parser.add_argument("--version", action="store_true", help="Display version information.")
	parser.add_argument("--manual", action="store_true", help="(Attempt to) open browser and show help")
	parser.add_argument("-name", "-n", help="'Stem' name for all output files.")
	parser.add_argument("-dna", "-d", help="Unaligned DNA (in FASTA format).")
	parser.add_argument("-gene", "-g", help="The genes to search for (multiple genes are comma-separated)")
	parser.add_argument("-nGenes", "-ng", help="The number of genes to search for (if fewer than suggested in 'genes' are required)")
	parser.add_argument("-species", "-s", help="Binomial names of species, each on a new line")
	parser.add_argument("-alignment", "-a", help="Aligned DNA (in FASTA format).")
	parser.add_argument("-constraint", "-c", help="Constraint tree (in newick format).")
	parser.add_argument("-email", "-e", help="Email address for GenBank searches.")
	parser.add_argument("-parameters", "-p", help="Parameter file giving detailed instructions to phyloGen.")
	main()
