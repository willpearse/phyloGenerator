#!/usr/bin/env python
# encoding: utf-8
"""
phyloGenerator.py
Created by Will Pearse on 2011-08-24.
Copyright (c) 2011 Imperial College london. All rights reserved.
TO-DO:
* BEAST
* 'Execute later' code
* Be able to delete a particular gene from the line-up
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
import pdb
def taxonIDLookup(taxonID):
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
	handleSearch = Entrez.esearch(db="taxonomy", term=spName)
	resultsSearch = Entrez.read(handleSearch)
	handleSearch.close()
	if resultsSearch['IdList']:
		return taxonIDLookup(resultsSearch['IdList'][0])
	else:
		return(tuple())

def cladeSpecies(cladeName):
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
	handle = Entrez.esearch(db="nucleotide",term=term, usehistory=usehistory, retStart=retStart, retMax=retMax)
	results = Entrez.read(handle)
	handle.close()
	return results

def eFetchSeqID(seqID, rettype='gb'):
	handle = Entrez.efetch(db="nucleotide", rettype=rettype, id=seqID)
	results = SeqIO.read(handle,rettype)
	handle.close()
	return results

def eFetchESearch(eSearchOutput, rettype='gb'):
	handle = Entrez.efetch(db="nucleotide", rettype=rettype, webenv=eSearchOutput['WebEnv'], query_key=eSearchOutput['QueryKey'])
	results = SeqIO.read(handle, rettype)
	handle.close()
	return results

def eSummary(seqID):
	handle = Entrez.esummary(db="nucleotide", id=seqID)
	results = Entrez.read(handle)
	handle.close()
	return results[0]

def sequenceDownload(spName, geneName=None, thorough=False, rettype='gb', titleText=None, noSeqs=1, seqChoice='random', download=True, retStart=0, retMax=20, targetLength=None, trimSeq=False, DNAtype='Standard', gapType='-', includeGenome=True, includePartial=True):
	def dwnSeq(includeGenome, includePartial):
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

def findGenes(speciesList, geneNames, download=False, titleText=None, targetNoGenes=-1, noSeqs=1, includePartial=True, includeGenome=True, seqChoice='random', verbose=True, thorough=False, spacer=10, delay=5):
	def findBestGene(foundSeqsArray):
		geneHits = foundSeqsArray.sum(axis=0)
		for i in range(len(geneHits)):
			if geneHits[i] == max(geneHits):
				return i
	
	if type(geneNames) is list:
		if targetNoGenes == -1:
			targetNoGenes = None
		
		#Download number of genes and histories for each species
		searchResults = []
		foundSeqs = []
		foundSeqsBool = np.zeros((len(speciesList), len(geneNames)), int)
		counter = 0
		for i in range(len(speciesList)):
			if verbose: print "Searching for:", speciesList[i]
			speciesGenes = []
			for k in range(len(geneNames)):
				sequence = sequenceDownload(speciesList[i], geneNames[k], titleText=titleText, noSeqs=noSeqs, includePartial=includePartial, includeGenome=includeGenome, seqChoice=seqChoice, download=download, thorough=thorough)
				counter += 1
				if counter == spacer:
					time.sleep(delay)
					counter = 0
				if download:
					speciesGenes.append(sequence)
				if sequence:
					foundSeqsBool[i,k] = 1
			if download:
				foundSeqs.append(speciesGenes)
		if targetNoGenes:
			currentFoundSeqs = foundSeqsBool[:]
			currentGeneNames = geneNames[:]
			bestGenes = []
			for each in range(targetNoGenes):
				bestGeneIndex = findBestGene(currentFoundSeqs)
				bestGenes.append(currentGeneNames[bestGeneIndex])
				currentFoundSeqs = np.delete(currentFoundSeqs, bestGeneIndex, axis=1)
				del currentGeneNames[bestGeneIndex]
			if download:
				output = []
				for i,sp in enumerate(speciesList):
					spOutput = []
					for j,gene in enumerate(geneNames):
						if gene in bestGenes:
							spOutput.append(foundSeqs[i][j])
					output.append(spOutput)
				return (output, bestGenes)
			else:
				return (currentFoundSeqs, bestGenes)
		else:
			if download:
				return (foundSeqs, geneNames)
			else:
				return (foundSeqsBool, geneNames)
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
	if method == 'quick': method = 'muscle-mafft-clustalo'
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
	
	def alignLen(gene):
		output = []
		for method in gene:
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

def RAxML(alignment, method='localVersion', tempStem='temp', outgroup=None, timeout=999999999, cladeList=None,  DNAmodel='GTR+G', constraint=None, cleanup=True, runNow=True):
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
	AlignIO.write(alignment, inputFile, "phylip-relaxed")
	#Constraint
	if constraint:
		if constraint.is_bifurcating():
			options += " -r " + tempStem + "_constraint.tre"
		else:
			options += " -g " + tempStem + "_constraint.tre"
		Phylo.write(constraint, tempStem + "_constraint.tre", 'newick')
		if constraint.total_branch_length() > 0:
			output = ''
			with open(tempStem + "_constraint.tre") as f:
				for each in f:
					output += each.strip()
			output = re.sub('[0-9]', '', output)
			output = re.sub('\.', '', output)
			output = re.sub('\:', '', output)
			with open(tempStem + "_constraint.tre", 'w') as f:
				f.write(output)
	commandLine = raxmlVersion + raxmlCompile + fileLine + algorithm + DNAmodel + options
	
	if not runNow:
		return commandLine
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
				if tempStem+"In.phylip.reduced"==each:
					os.remove(each)
		return tree
			
	else:
		raise RuntimeError("Either phylogeny building program failed, or ran out of time")

def BEAST(alignment, method='GTR+GAMMA', tempStem='temp', timeout=999999999, constraint=None, cleanup=True, runNow=True, chainLength=1000000, logRate=1000, screenRate=1000):
	completeConstraint = False
	with open(tempStem+"_BEAST.xml", 'w') as f:
		f.write('<?xml version="1.0" standalone="yes"?>\n')
		f.write('<beast>\n')
		f.write('	<!-- The list of taxa analyse (can also include dates/ages).				 -->\n')
		f.write('	<taxa id="taxa">\n')
		for each in alignment:
			f.write('		<taxon id="'+each.id+'"/>\n')
		f.write('	</taxa>\n')
		if constraint:
			clades = []
			for clade in constraint.find_clades():
				temp = [x.name for x in clade.get_terminals()]
				if len(temp) > 1:
					clades.append(temp)
			if len(set([item for sublist in clades for item in sublist])) == len(constraint.get_terminals()):
				completeConstraint = True
			for i,clade in enumerate(clades):
				f.write('	<taxa id="cladeNo' + str(i) + '">\n')
				for sp in clade:
					f.write('		<taxon idref="' + sp + '"/>\n')
				f.write('	</taxa>')
		f.write('	<!-- The sequence alignment (each sequence refers to a taxon above).		 -->\n')
		f.write('	<alignment id="alignment" dataType="nucleotide">\n')
		for each in alignment:
			f.write('		<sequence>\n')
			f.write('			<taxon idref="'+each.id+'"/>\n')
			f.write(each.seq.tostring()+'\n')
			f.write('		</sequence>\n')
		f.write('	</alignment>\n')
		f.write('	<!-- The unique patterns from 1 to end										 -->\n')
		f.write('	<patterns id="patterns" from="1">\n')
		f.write('		<alignment idref="alignment"/>\n')
		f.write('	</patterns>\n')
		f.write('	<!-- A prior on the distribution node heights defined given					 -->\n')
		f.write('	<!-- a Yule speciation process (a pure birth process).						 -->\n')
		f.write('	<yuleModel id="yule" units="substitutions">\n')
		f.write('		<birthRate>\n')
		f.write('			<parameter id="yule.birthRate" value="1.0" lower="0.0" upper="Infinity"/>\n')
		f.write('		</birthRate>\n')
		f.write('	</yuleModel>\n')
		f.write('	<!-- This is a simple constant population size coalescent model				 -->\n')
		f.write('	<!-- that is used to generate an initial tree for the chain.				 -->\n')
		f.write('	<constantSize id="initialDemo" units="substitutions">\n')
		f.write('		<populationSize>\n')
		f.write('			<parameter id="initialDemo.popSize" value="100.0"/>\n')
		f.write('		</populationSize>\n')
		f.write('	</constantSize>\n')
		if completeConstraint:
			f.write('	<newick id="startingTree">')
			Phylo.write(constraint, 'tempStem'+'_TREE.tre', 'newick')
			with open('tempStem'+'_TREE.tre') as tF:
				treeFormat = tF.readlines()[0]
			os.remove()
			f.write('		'+treeFormat)
			f.write('	</newick>')
		else:
			f.write('	<!-- Generate a random starting tree under the coalescent process			 -->\n')
			f.write('	<coalescentTree id="startingTree" rootHeight="0.092">\n')
			f.write('		<taxa idref="taxa"/>\n')
			f.write('		<constantSize idref="initialDemo"/>\n')
			f.write('	</coalescentTree>\n')
		f.write('	<!-- Generate a tree model													 -->\n')
		f.write('	<treeModel id="treeModel">\n')
		if completeConstraint:
			f.write('		<newick idref="startingTree"/>')
		else:
			f.write('		<coalescentTree idref="startingTree"/>\n')
		f.write('		<rootHeight>\n')
		f.write('			<parameter id="treeModel.rootHeight"/>\n')
		f.write('		</rootHeight>\n')
		f.write('		<nodeHeights internalNodes="true">\n')
		f.write('			<parameter id="treeModel.internalNodeHeights"/>\n')
		f.write('		</nodeHeights>\n')
		f.write('		<nodeHeights internalNodes="true" rootNode="true">\n')
		f.write('			<parameter id="treeModel.allInternalNodeHeights"/>\n')
		f.write('		</nodeHeights>\n')
		f.write('	</treeModel>\n')
		f.write('	<!-- Generate a speciation likelihood for Yule or Birth Death				 -->\n')
		f.write('	<speciationLikelihood id="speciation">\n')
		f.write('		<model>\n')
		f.write('			<yuleModel idref="yule"/>\n')
		f.write('		</model>\n')
		f.write('		<speciesTree>\n')
		f.write('			<treeModel idref="treeModel"/>\n')
		f.write('		</speciesTree>\n')
		f.write('	</speciationLikelihood>\n')
		f.write('	<!-- The uncorrelated relaxed clock (Drummond, Ho, Phillips & Rambaut, 2006) -->\n')
		f.write('	<discretizedBranchRates id="branchRates">\n')
		f.write('		<treeModel idref="treeModel"/>\n')
		f.write('		<distribution>\n')
		f.write('			<logNormalDistributionModel meanInRealSpace="true">\n')
		f.write('				<mean>\n')
		f.write('					<parameter id="ucld.mean" value="1.0"/>\n')
		f.write('				</mean>\n')
		f.write('				<stdev>\n')
		f.write('					<parameter id="ucld.stdev" value="0.3333333333333333" lower="0.0" upper="Infinity"/>\n')
		f.write('				</stdev>\n')
		f.write('			</logNormalDistributionModel>\n')
		f.write('		</distribution>\n')
		f.write('		<rateCategories>\n')
		f.write('			<parameter id="branchRates.categories" dimension="18"/>\n')
		f.write('		</rateCategories>\n')
		f.write('	</discretizedBranchRates>\n')
		f.write('	<rateStatistic id="meanRate" name="meanRate" mode="mean" internal="true" external="true">\n')
		f.write('		<treeModel idref="treeModel"/>\n')
		f.write('		<discretizedBranchRates idref="branchRates"/>\n')
		f.write('	</rateStatistic>\n')
		f.write('	<rateStatistic id="coefficientOfVariation" name="coefficientOfVariation" mode="coefficientOfVariation" internal="true" external="true">\n')
		f.write('		<treeModel idref="treeModel"/>\n')
		f.write('		<discretizedBranchRates idref="branchRates"/>\n')
		f.write('	</rateStatistic>\n')
		f.write('	<rateCovarianceStatistic id="covariance" name="covariance">\n')
		f.write('		<treeModel idref="treeModel"/>\n')
		f.write('		<discretizedBranchRates idref="branchRates"/>\n')
		f.write('	</rateCovarianceStatistic>\n')
		if 'GTR' in method:
			f.write('	<!-- The general time reversible (GTR) substitution model					 -->\n')
			f.write('	<gtrModel id="gtr">\n')
			f.write('		<frequencies>\n')
			f.write('			<frequencyModel dataType="nucleotide">\n')
			f.write('				<frequencies>\n')
			if 'base=empirical' in method:
				f.write('					<parameter id="frequencies" dimension="4"/>\n')
			else:
				f.write('					<parameter id="frequencies" value="0.25 0.25 0.25 0.25"/>\n')
			f.write('				</frequencies>\n')
			f.write('			</frequencyModel>\n')
			f.write('		</frequencies>\n')
			f.write('		<rateAC>\n')
			f.write('			<parameter id="ac" value="1.0" lower="0.0" upper="Infinity"/>\n')
			f.write('		</rateAC>\n')
			f.write('		<rateAG>\n')
			f.write('			<parameter id="ag" value="1.0" lower="0.0" upper="Infinity"/>\n')
			f.write('		</rateAG>\n')
			f.write('		<rateAT>\n')
			f.write('			<parameter id="at" value="1.0" lower="0.0" upper="Infinity"/>\n')
			f.write('		</rateAT>\n')
			f.write('		<rateCG>\n')
			f.write('			<parameter id="cg" value="1.0" lower="0.0" upper="Infinity"/>\n')
			f.write('		</rateCG>\n')
			f.write('		<rateGT>\n')
			f.write('			<parameter id="gt" value="1.0" lower="0.0" upper="Infinity"/>\n')
			f.write('		</rateGT>\n')
			f.write('	</gtrModel>\n')
		elif 'HKY' in method:
			f.write('	<!-- The HKY substitution model (Hasegawa, Kishino & Yano, 1985)             -->')
			f.write('	<HKYModel id="hky">')
			f.write('		<frequencies>')
			f.write('			<frequencyModel dataType="nucleotide">')
			f.write('				<frequencies>')
			f.write('					<parameter id="frequencies" value="0.25 0.25 0.25 0.25"/>')
			f.write('				</frequencies>')
			f.write('			</frequencyModel>')
			f.write('		</frequencies>')
			f.write('		<kappa>')
			f.write('			<parameter id="kappa" value="2.0" lower="0.0" upper="Infinity"/>')
			f.write('		</kappa>')
			f.write('	</HKYModel>')
		else:
			raise RuntimeError("No valid DNA substituion model specified for BEAST.")
		f.write('	<!-- site model																 -->\n')
		f.write('	<siteModel id="siteModel">\n')
		f.write('		<substitutionModel>\n')
		if 'GTR' in method:
			f.write('			<gtrModel idref="gtr"/>\n')
		elif 'HKY' in method:
			f.write('			<HKYModel idref="hky"/>')
		else:
			raise RuntimeError("No valid DNA substituion model specified for BEAST.")
		f.write('		</substitutionModel>\n')
		if 'GAMMA' in method:
			f.write('		<gammaShape gammaCategories="4">')
			f.write('			<parameter id="alpha" value="0.5" lower="0.0" upper="1000.0"/>')
			f.write('		</gammaShape>')
		f.write('	</siteModel>\n')
		f.write('	<treeLikelihood id="treeLikelihood" useAmbiguities="false">\n')
		f.write('		<patterns idref="patterns"/>\n')
		f.write('		<treeModel idref="treeModel"/>\n')
		f.write('		<siteModel idref="siteModel"/>\n')
		f.write('		<discretizedBranchRates idref="branchRates"/>\n')
		f.write('	</treeLikelihood>\n')
		f.write('	<!-- Define operators														 -->\n')
		f.write('	<operators id="operators">\n')
		if 'GTR' in method:
			f.write('		<scaleOperator scaleFactor="0.75" weight="0.1">\n')
			f.write('			<parameter idref="ac"/>\n')
			f.write('		</scaleOperator>\n')
			f.write('		<scaleOperator scaleFactor="0.75" weight="0.1">\n')
			f.write('			<parameter idref="ag"/>\n')
			f.write('		</scaleOperator>\n')
			f.write('		<scaleOperator scaleFactor="0.75" weight="0.1">\n')
			f.write('			<parameter idref="at"/>\n')
			f.write('		</scaleOperator>\n')
			f.write('		<scaleOperator scaleFactor="0.75" weight="0.1">\n')
			f.write('			<parameter idref="cg"/>\n')
			f.write('		</scaleOperator>\n')
			f.write('		<scaleOperator scaleFactor="0.75" weight="0.1">\n')
			f.write('			<parameter idref="gt"/>\n')
			f.write('		</scaleOperator>\n')
		elif 'HKY' in method:
			f.write('		<deltaExchange delta="0.01" weight="0.1">\n')
			f.write('			<parameter idref="frequencies"/>\n')
			f.write('		</deltaExchange>\n')
		else:
			raise RuntimeError("No valid DNA substituion model specified for BEAST.")
		if 'GAMMA' in method:
			f.write('		<scaleOperator scaleFactor="0.75" weight="0.1">')
			f.write('			<parameter idref="alpha"/>')
			f.write('		</scaleOperator>')
		f.write('		<scaleOperator scaleFactor="0.75" weight="3">\n')
		f.write('			<parameter idref="ucld.stdev"/>\n')
		f.write('		</scaleOperator>\n')
		if not completeConstraint:
			f.write('		<subtreeSlide size="0.0092" gaussian="true" weight="15">\n')
			f.write('			<treeModel idref="treeModel"/>\n')
			f.write('		</subtreeSlide>\n')
			f.write('		<narrowExchange weight="15">\n')
			f.write('			<treeModel idref="treeModel"/>\n')
			f.write('		</narrowExchange>\n')
			f.write('		<wideExchange weight="3">\n')
			f.write('			<treeModel idref="treeModel"/>\n')
			f.write('		</wideExchange>\n')
			f.write('		<wilsonBalding weight="3">\n')
			f.write('			<treeModel idref="treeModel"/>\n')
			f.write('		</wilsonBalding>\n')
			f.write('		<scaleOperator scaleFactor="0.75" weight="3">\n')
			f.write('			<parameter idref="treeModel.rootHeight"/>\n')
			f.write('		</scaleOperator>\n')
		f.write('		<uniformOperator weight="30">\n')
		f.write('			<parameter idref="treeModel.internalNodeHeights"/>\n')
		f.write('		</uniformOperator>\n')
		f.write('		<scaleOperator scaleFactor="0.75" weight="3">\n')
		f.write('			<parameter idref="yule.birthRate"/>\n')
		f.write('		</scaleOperator>\n')
		f.write('		<upDownOperator scaleFactor="0.75" weight="3">\n')
		f.write('			<up>\n')
		f.write('			</up>\n')
		f.write('			<down>\n')
		f.write('				<parameter idref="treeModel.allInternalNodeHeights"/>\n')
		f.write('			</down>\n')
		f.write('		</upDownOperator>\n')
		f.write('		<swapOperator size="1" weight="10" autoOptimize="false">\n')
		f.write('			<parameter idref="branchRates.categories"/>\n')
		f.write('		</swapOperator>\n')
		f.write('		<randomWalkIntegerOperator windowSize="1" weight="10">\n')
		f.write('			<parameter idref="branchRates.categories"/>\n')
		f.write('		</randomWalkIntegerOperator>\n')
		f.write('		<uniformIntegerOperator weight="10">\n')
		f.write('			<parameter idref="branchRates.categories"/>\n')
		f.write('		</uniformIntegerOperator>\n')
		f.write('	</operators>\n')
		f.write('	<!-- Define MCMC															 -->\n')
		f.write('	<mcmc id="mcmc" chainLength="' + str(chainLength) + '" autoOptimize="true">\n')
		f.write('		<posterior id="posterior">\n')
		f.write('			<prior id="prior">\n')
		if 'GTR' in method:
			f.write('				<gammaPrior shape="0.05" scale="10.0" offset="0.0">\n')
			f.write('					<parameter idref="ac"/>\n')
			f.write('				</gammaPrior>\n')
			f.write('				<gammaPrior shape="0.05" scale="20.0" offset="0.0">\n')
			f.write('					<parameter idref="ag"/>\n')
			f.write('				</gammaPrior>\n')
			f.write('				<gammaPrior shape="0.05" scale="10.0" offset="0.0">\n')
			f.write('					<parameter idref="at"/>\n')
			f.write('				</gammaPrior>\n')
			f.write('				<gammaPrior shape="0.05" scale="10.0" offset="0.0">\n')
			f.write('					<parameter idref="cg"/>\n')
			f.write('				</gammaPrior>\n')
			f.write('				<gammaPrior shape="0.05" scale="10.0" offset="0.0">\n')
			f.write('					<parameter idref="gt"/>\n')
			f.write('				</gammaPrior>\n')
		f.write('				<exponentialPrior mean="0.3333333333333333" offset="0.0">\n')
		f.write('					<parameter idref="ucld.stdev"/>\n')
		f.write('				</exponentialPrior>\n')
		f.write('				<speciationLikelihood idref="speciation"/>\n')
		f.write('			</prior>\n')
		f.write('			<likelihood id="likelihood">\n')
		f.write('				<treeLikelihood idref="treeLikelihood"/>\n')
		f.write('			</likelihood>\n')
		f.write('		</posterior>\n')
		f.write('		<operators idref="operators"/>\n')
		f.write('		<!-- write log to screen													 -->\n')
		f.write('		<log id="screenLog" logEvery="' + str(screenRate) + '">\n')
		f.write('			<column label="Posterior" dp="4" width="12">\n')
		f.write('				<posterior idref="posterior"/>\n')
		f.write('			</column>\n')
		f.write('			<column label="Prior" dp="4" width="12">\n')
		f.write('				<prior idref="prior"/>\n')
		f.write('			</column>\n')
		f.write('			<column label="Likelihood" dp="4" width="12">\n')
		f.write('				<likelihood idref="likelihood"/>\n')
		f.write('			</column>\n')
		f.write('			<column label="rootHeight" sf="6" width="12">\n')
		f.write('				<parameter idref="treeModel.rootHeight"/>\n')
		f.write('			</column>\n')
		f.write('			<column label="ucld.mean" sf="6" width="12">\n')
		f.write('				<parameter idref="ucld.mean"/>\n')
		f.write('			</column>\n')
		f.write('		</log>\n')
		f.write('		<!-- write log to file														 -->\n')
		f.write('		<log id="fileLog" logEvery="' + str(logRate) + '" fileName="'+tempStem+'.log" overwrite="false">\n')
		f.write('			<posterior idref="posterior"/>\n')
		f.write('			<prior idref="prior"/>\n')
		f.write('			<likelihood idref="likelihood"/>\n')
		f.write('			<parameter idref="treeModel.rootHeight"/>\n')
		f.write('			<parameter idref="yule.birthRate"/>\n')
		if 'GTR' in method:
			f.write('			<parameter idref="ac"/>\n')
			f.write('			<parameter idref="ag"/>\n')
			f.write('			<parameter idref="at"/>\n')
			f.write('			<parameter idref="cg"/>\n')
			f.write('			<parameter idref="gt"/>\n')
			f.write('			<parameter idref="frequencies"/>\n')
		elif 'HKY' in method:
			f.write('			<parameter idref="kappa"/>')
		else:
			raise RuntimeError("No valid DNA substituion model specified for BEAST.")
		if 'GAMMA' in method:
			f.write('			<parameter idref="alpha"/>')
		f.write('			<parameter idref="ucld.mean"/>\n')
		f.write('			<parameter idref="ucld.stdev"/>\n')
		f.write('			<rateStatistic idref="meanRate"/>\n')
		f.write('			<rateStatistic idref="coefficientOfVariation"/>\n')
		f.write('			<rateCovarianceStatistic idref="covariance"/>\n')
		f.write('			<treeLikelihood idref="treeLikelihood"/>\n')
		f.write('			<speciationLikelihood idref="speciation"/>\n')
		f.write('		</log>\n')
		f.write('		<!-- write tree log to file													 -->\n')
		f.write('		<logTree id="treeFileLog" logEvery="' + str(logRate) + '" nexusFormat="true" fileName="'+tempStem+'.trees" sortTranslationTable="true">\n')
		f.write('			<treeModel idref="treeModel"/>\n')
		f.write('			<discretizedBranchRates idref="branchRates"/>\n')
		f.write('			<posterior idref="posterior"/>\n')
		f.write('		</logTree>\n')
		f.write('	</mcmc>\n')
		f.write('	<report>\n')
		f.write('		<property name="timer">\n')
		f.write('			<mcmc idref="mcmc"/>\n')
		f.write('		</property>\n')
		f.write('	</report>\n')
		f.write('</beast>\n')
	
	commandLine = "beast " + tempStem + "_BEAST.xml"
	
	if not runNow:
		return commandLine
	
	pipe = TerminationPipe(commandLine, timeout)
	pipe.run(silent=False)
	pdb.set_trace()
	if not pipe.failure:
		print "...removing burn-in of 10%..."
		commandLine = 'treeannotator -burnin 1000 -height median ' + tempStem + ".trees" + tempStem + "Final.tre"
		pipeAnotate = TerminationPipe(commandLine, timeout)
		pipeAnotate.run()
		if not pipeAnotate.failure:
			tree = Phylo.read(tempStem + "Final.tre", "newick")
			if cleanup:
				os.remove(tempStem + ".trees")
				os.remove(tempStem + ".log")
			return tree
		else:
			raise RuntimeError("Either tree annotation failed, or ran out of time")
	else:
		raise RuntimeError("Either phylogeny building program failed, or ran out of time")

def trimSequence(seq, DNAtype='Standard', gapType='-'):
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
		taxaFile = ' -t ' + spNames
		commandLine = 'phylomatic' + phylogenyFile + taxaFile + " > phylomaticPhyloGenerator" + tempStem
		pipe = TerminationPipe(commandLine, 999)
		pipe.run(silent=True)
		if not pipe.failure:
			constraint = Phylo.read("phylomaticPhyloGenerator" + tempStem, 'newick')
			os.remove("phylomaticPhyloGenerator" + tempStem)
			return constraint
		else:
			raise RuntimeError("Phylomatic did not run correctly")
		
	elif method == "GenBank":
		print "THIS DOESN'T WORK!!!!!!!!"
		lineages = [findLineage(x) for x in spNames]
		return recursiveTree(lineages)
	else:
		raise RuntimeError("Unrecognised constraint tree creation method specified")

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
			tStdout	 = open('termPipeStdOut.txt', 'w')
			tStderr	 = open('termPipeStdErr.txt', 'w')
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
	def __init__(self, args):
		self.fastaFile = False
		self.GenBankFile = False
		self.sequences = []
		self.speciesNames = []
		self.dnaCheck = []
		self.phylogeny = False
		self.alignment = []
		self.smoothPhylogeny = False
		self.root = False
		self.genBankIDs = []
		self.constraint = False
		self.genBankRetMax = 50
		self.spacer = 10
		self.alignmentMethod = []
		self.alignmentMethodChosen = []
		self.email = ''
		self.codonModels = []
		self.taxonomy = []
		self.raxml = 'RAxML-localVersion'
		self.initialSeqChoice = 'medianLength'
		self.replaceSeqChoice = 'targetLength'
		self.uniqueTaxonomy = []
		self.tracker = 0
		
		#Stem name
		if args.name:
			self.stem = args.name
			print "\nUsing stem name", args.name, "..."
		else:
			print "\nPlease input a 'stem' name for all your output (phylogeny, sequences, etc.)"
			self.stem = raw_input("Stem name: ")
		
		#Email
		if args.email:
			self.email = args.email
			Entrez.email = self.email
		
		#Timeout
		if args.delay:
			self.delay = int(args.delay)
		else:
			self.delay = 5
		
		#Gene
		self.nGenes = -1
		if args.gene:
			self.genes = args.gene.split(',')
			print "Using gene(s)", self.genes
		else:
			print"Please enter the name of the gene you want to use, e.g. 'COI' for cytochrome oxidase one. To enter multiple genes, enter each on a separateline."
			locker = True
			self.genes = []
			while locker:
				inputGene = raw_input("")
				if inputGene:
					self.genes.append(inputGene)
					self.codonModels.append('Standard')
				else:
					if self.genes:
						locker = False
					else:
						print "...you cannot continue without enter at least one gene name."
		
		#Gene Number
		if args.nGenes:
			self.nGenes = int(args.nGenes)
		
		#FASTA input
		if args.dna:
			self.loadDNAFile(args.dna)
		
		if args.species:
			self.loadGenBank(args.species)
		
		#Alignment method
		# - stupidly, this contrasts with 'self.alignmentMethods', which is a list... Fix this when you clean up the self.align method...
		if args.alignment:
			self.alignmentMethod = args.alignment
		
		#Constraint tree
		self.constraintMethod = ''
		#Phylomatic
		if args.consPhylomat:
			self.phylomaticPhylogeny, self.phylomaticTaxonomy = args.consPhylom.split(",")
			self.constraintMethod = 'phylomatic'
		
		#Pre-supplied
		if args.consTree:
			self.constraintFilename = args.consTree
			self.constraintMethod = 'newick'
		
		#Options file
		if args.parameters:
			options = []
			try:
				with open(args.parameters) as f:
					for each in f:
						options.append(each.strip())
			except:
				print "\nERROR: Could not find parameters file. Exiting."
			for line in options:
				if 'GenBank.retmax' in line:
					self.genBankRetMax = int(line.replace("GenBank.retmax=", ""))
				elif 'GenBank.spacer' in line:
					self.spacer = int(line.replace("GenBank.spacer=", ""))
				elif 'RAxML' in line:
					self.raxml = line
				elif 'GenBank.initialSeqChoice' in line:
					self.initialSeqChoice = line.replace("GenBank.initialSeqChoice=", "")
				elif 'GenBank.replaceSeqChoice' in line:
					temp = line.replace("GenBank.repalceSeqChoice=", "")
					if temp == 'geneMedianLength':
						self.replaceSeqChoice = 'targetLength'
					else:
						self.replaceSeqChoice = temp
	
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
				self.sequences, self.genes = findGenes(self.speciesNames, self.genes, seqChoice=self.initialSeqChoice, verbose=True, download=True, thorough=True, targetNoGenes=self.nGenes, spacer=self.spacer, delay=self.delay)
	
	def dnaChecking(self, tolerance=0.1):
		self.tolerance = tolerance
		self.dnaCheck = checkSequenceList(self.sequences, tolerance=self.tolerance, method="quantileDetect")
		sequenceDisplay(self.sequences, self.speciesNames, self.genes, self.dnaCheck)
	
	def dnaEditing(self):
		def deleteMode(firstTime=True):
			if firstTime:
				print "\nYou're in deletion mode. To delete a species, enter its SeqID and press return."
				print "\tTo delete an entire gene, type 'gene', press enter, then the name of the gene you want to delete."
				print "\t*One species at a time please!*"
				print "To change to the 'reload', 'trim' or 'replace' modes, simply type their names then hit enter."
			inputSeq = raw_input("DNA Editing (delete): ")
			if inputSeq:
				try:
					if int(inputSeq) in range(len(self.speciesNames)):
						del self.sequences[int(inputSeq)]
						del self.speciesNames[int(inputSeq)]
						print "SeqID", inputSeq, "Successfully deleted"
						print "Re-calulating summary statistics..."
						self.dnaChecking()
						return 'delete', False
				except:
					pass
				if inputSeq == "gene":
					geneLocker = True
					while geneLocker:
						print "Enter the name of the gene you want to delete, or just hit enter to cancel."
						geneInput = raw_input("DNA Editing (delete gene): ")
						if geneInput:
							if geneInput in self.genes:
								for i,species in enumerate(self.sequences):
									for j,gene in enumerate(species):
										if self.genes[j] == geneInput:
											del self.sequences[i][j]
											continue
								for i,gene in enumerate(self.genes):
									if gene == geneInput:
										del self.genes[i]
										break
								if self.nGenes != -1:
									self.nGenes -= 1
								geneLocker = False
								print "Gene", geneInput, "successfully deleted. Re-calculating summary statistics..."
								self.dnaChecking()
								return 'delete', False
							else:
								print "Sorry,", inputSeq, "was not recognised. Please try again."
						else:
							return 'delete', False
				elif inputSeq == "trim":
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
					if int(inputSeq) in range(len(self.speciesNames)):
						for i,each in enumerate(self.sequences[int(inputSeq)]):
							self.APICheck()
							self.sequences[int(inputSeq)][i] = sequenceDownload(self.speciesNames[int(inputSeq)], self.genes[i], thorough=True, retMax=self.maxGenBankDownload, seqChoice=self.replaceSeqChoice, targetLength=self.dnaCheck['quantileLengths'][i][2])
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
							self.APICheck()
							self.sequences[seqID][i] = sequenceDownload(self.speciesNames[seqID], self.genes[i], thorough=True, retMax=self.maxGenBankDownload, seqChoice=self.replaceSeqChoice, targetLength=self.dnaCheck['quantileLengths'][i][2])
							print "Re-calulating summary statistics..."
							self.dnaChecking()
							return 'reload', False
				if ">" in inputSeq:
					threshold = int(inputSeq.replace('>', ''))
					for i,sp in enumerate(self.sequences):
						for j,gene in enumerate(sp):
							if len(self.sequences[i][j]) > threshold:
								self.APICheck()
								self.sequences[i][j] = sequenceDownload(self.speciesNames[i], self.genes[j], thorough=True, retMax=self.maxGenBankDownload, seqChoice=self.replaceSeqChoice, targetLength=self.dnaCheck['quantileLengths'][j][2])
					print "Re-calulating summary statistics..."
					self.dnaChecking()
					return 'reload', False
				if "<" in inputSeq:
					threshold = int(inputSeq.replace('<', ''))
					for i,sp in enumerate(self.sequences):
						for j,gene in enumerate(sp):
							if len(self.sequences[i][j]) < threshold:
								self.APICheck()
								self.sequences[i][j] = sequenceDownload(self.speciesNames[j], self.genes[i], thorough=True, retMax=self.maxGenBankDownload, seqChoice=self.replaceSChoice, targetLength=self.dnaCheck['quantileLengths'][j][2])
					print "Re-calulating summary statistics..."
					self.dnaChecking()
					return 'reload', False
				if inputSeq == "EVERYTHING":
					geneOutput = findGenes(self.speciesNames, self.genes, seqChoice=self.replaceSeqChoice, verbose=True, thorough=True, download=True, retMax=self.maxGenBankDownload, targetNoGenes=self.nGenes, delay=self.delay, spacer=self.spacer)
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
					if int(inputSeq) in range(len(self.speciesNames)):
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
				print "\tTo replace all species without any sequences with a relative it is more related to than any other species according to the GenBank taxonomy, type 'THOROUGH' and press enter."
				print "To change to the 'delete', 'trim' or 'reload' modes, simply type their names then hit enter."
			inputSeq = raw_input("DNA Editing (replace):")
			if inputSeq:
				try:
					if int(inputSeq) in range(len(self.speciesNames)):
						i = int(inputSeq)
						lineage = findLineage(self.speciesNames[i])
						if lineage:
							replacements = cladeSpecies(lineage[1])
							print "...looking for alternatives for", self.speciesNames[i], "in clade", lineage[1]
							for candidate in replacements:
								temp = []
								locker = False
								for j,gene in enumerate(self.genes):
									self.APICheck()
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
								self.APICheck()
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
					for i,sp in enumerate(self.speciesNames):
						tracker = 0
						for j,gene in enumerate(sp):
							if gene:
								tracker += 1
						if tracker == 0:
							lineage = findLineage(self.speciesNames[i])
							self.APICheck()
							replacements = cladeSpecies(lineage[1])
							print "...looking for alternatives for", self.speciesNames[i], "in clade", lineage[1]
							for candidate in replacements:
								temp = []
								locker = False
								for gene in self.genes:
									self.APICheck()
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
				elif inputSeq == "THOROUGH":
					#Find taxonomy
					for i,sp in enumerate(self.speciesNames):
						self.taxonomy.append(findLineage(sp))
					
					print "...Taxonomy downloaded"
					#Find highest unique level for each species
					for i,sp in enumerate(self.speciesNames):
						uniqueLineage = []
						finished = False
						for currLevel in self.taxonomy[i]:
							for j,spp in enumerate(self.speciesNames):
								if i!=j:
									if currLevel in self.taxonomy[j]:
										break
							else:
								uniqueLineage.append(currLevel)
						self.uniqueTaxonomy.append(uniqueLineage)
					
					print "...Taxonomy unique to each species found\n"
					
					#Find any replacements, if necessary, using the taxonomy
					for i,sp in enumerate(self.speciesNames):
						tracker = 0
						for j,gene in enumerate(self.genes):
							if self.sequences[i][j]:
								tracker += 1
						if tracker == 0:
							print "...looking for alternative for", sp
							if len(self.uniqueTaxonomy[i]) > 1:
								pdb.set_trace()
								for candidate in self.uniqueTaxonomy[i][1:]:
									geneList = []
									locker = False
									for gene in self.genes:
										self.APICheck()
										temp = sequenceDownload(candidate, gene)
										geneList.append(temp)
										if temp:
											self.speciesNames[i] = candidate[0] +"(" + self.speciesNames[i] + ")"
											locker = True
									if locker:
										self.sequences[i] = geneList
										print "......alternative found:", candidate
										break
							print "......can't find one"
					print "Re-calulating summary statistics..."
					self.dnaChecking()
					return 'replace', False
				elif inputSeq == "delete":
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
		alignmentDisplay(self.alignment, self.alignmentMethods, self.genes, checkAlignmentList(self.alignment, method='everything'))
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
				try:
					if int(choice) in range(len(self.alignmentMethods)):
						self.alignment[i] = self.alignment[i][int(choice)]
						self.alignmentMethodChosen.append(self.alignmentMethod[int(choice)])
						locker = False
					else:
						print "Sorry, didn't get that. Try again."
				except:
					print "Sorry, didn't get that. Try again."
	
	def align(self):
		methods = ['muscle', 'mafft', 'clustalo', 'prank', 'everything', 'quick']
		if not self.alignmentMethod:
			print "Enter the name of an alignment method ('muscle', 'mafft', 'clustalo', 'prank'), 'everything' to do all four and compare their outputs, 'quick' to do only the first three (prank can take a while!), or simply hit return to align with muscle."
			locker = True
			while locker:
				alignInput = raw_input("")
				if alignInput:
					if alignInput in methods:
						self.alignmentMethod = alignInput
						locker = False
					else:
						print "Sorry, I didn't recognise", alignInput, "- please try again."
				else:
					self.alignmentMethod = 'muscle'
					locker = False
		print "Starting alignment..."
		self.alignment = alignSequences(self.sequences, method=self.alignmentMethod, tempStem='temp', timeout=99999999, nGenes=len(self.genes))
		print "\nAlignment complete!"
		if self.alignmentMethod == 'everything':
			self.alignmentMethods = ['muscle', 'mafft', 'clustalo', 'prank']
		elif self.alignmentMethod == 'quick':
			self.alignmentMethods = ['muscle', 'mafft', 'clustalo']
		else:
			self.alignmentMethods = [self.alignmentMethod]
	
	def phylogen(self, method="RAxML-localVersion"):
		print"\n Running with options:", method
		self.phylogeny = RAxML(self.alignment, method=self.raxml, constraint=self.constraint, timeout=999)
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
			for k in range(len(self.genes)):
				if self.sequences[i][k]:
					tGenBankIDs.append(self.sequences[i][k].name)
					self.sequences[i][k].name = self.speciesNames[i].replace(" ", "_")
				else:
					tGenBankIDs.append("NO_SEQUENCE")
			self.genBankIDs.append(tGenBankIDs)
	
	def writeOutput(self):
		#Log - TO-DO
		#Sequences
		if self.sequences:
			for i,gene in enumerate(self.genes):
				currentGene = []
				for seq in self.sequences:
					currentGene.append(seq[i])
				SeqIO.write(currentGene, self.stem+"_"+gene+".fasta", 'fasta')
		
		#Alignment
		AlignIO.write(self.alignment, self.stem+"_alignment.fasta", 'fasta')
		
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
		
		#Taxonomy
		if self.taxonomy:
			with open(self.stem+"_taxonomy.txt", 'w') as f:
				f.write("Species Name, Taxonomy...\n")
				for i,each in enumerate(self.taxonomy):
					if each:
						f.write(self.speciesNames[i] + "," + ",".join(each) + "\n")
					else:
						f.write(self.speciesNames[i] + ", NOTHING FOUND\n")
		
		#Smoothed phylogeny
		if self.smoothPhylogeny:
			Phylo.write(self.smoothPhylogeny, self.stem+"_phylogeny_smoothed.tre", 'newick')
	
	def getConstraint(self, fileName=""):
		def newickConstraint():
			if self.constraintFile:
				try:
					self.constraint = Phylo.read(self.constraintFile, 'newick')
				except IOError:
					print "\nConstraint tree not found. Error!"
					sys.exit()
			else:
				print "\nEnter the filename of your constraint tree (in newick format) below. Just press enter to continue without a constraint tree."
				locker = True
				while locker:
					inputConstraint = raw_input("Newick Constraint Tree: ")
					if inputConstraint:
						try:
							self.constraint = Phylo.read(inputConstraint, 'newick')
						except IOError:
							print "\nFile not found. Please try again!"
						if self.checkConstraint():
							print "Constraint tree loaded!"
							locker = False
							return False
						else:
							self.constraint = False
							print "Constraint tree does *not* match the species names you've inputted. Please load another file."
					else:
						print "...No constraint tree loaded"
						locker = False
						return False
		
		def phylomatic():
			if self.phylomaticTaxonomy and self.phylomaticPhylogeny:
				try:
					self.constraint = createConstraintTree(self.phylomaticTaxonomy, method="phylomaticTaxonomy", fileName=self.phylomaticPhylogeny, tempStem='temp')
					return False
				except:
					print "Phylomatic constraint tree construction problem. Exiting."
					sys.exit()
			else:
				print "\nCreating a constraint tree using Phylomatic."
				if not self.phylomaticPhylogeny:
					print " Enter the filename of the reference phylogeny (in newick format) below. Just hit enter to carry on without a constraint tree."
					locker = True
					while locker:
						inputPhylogeny = raw_input("Phylomatic (reference tree): ")
						if inputPhylogeny:
							try:
								open(inputPhylogeny)
								print "...Reference tree found!"
								self.phylomaticPhylogeny = inputPhylogeny
								locker = False
							except IOError:
								print "\nFile not found. Please try again!"
						else:
							print "...No constraint tree loaded. Continuing."
							return False
				else:
					print "...Phlyomatic reference tree loaded!"
				if not self.phylomaticTaxonomy:
					print "Now enter the taxonomy (in Phylomatic's format) for the species you're using."
					locker = True
					while locker:
						inputTaxonomy = raw_input("Phylomatic (taxonomy): ")
						if inputTaxonomy:
							try:
								open(inputTaxonomy)
								print "...Reference taxonomy found!"
								self.phylomaticTaxonomy = inputTaxonomy
								locker = False
							except IOError:
								print "\nFile not found. Please try again!"
						else:
							print "...No constraint tree loaded. Continuing."
							locker = False
				else:
					print "...Phylomatic reference taxonomy loaded!"
				print "Attempting to run Phylomatic..."
				try:
					self.constraint = createConstraintTree(self.phylomaticTaxonomy, method="phylomaticTaxonomy", fileName=self.phylomaticPhylogeny, tempStem='temp')
					return False
				except:
					print "...something went wrong with that. Check your taxonomy and your phylogeny"
		
		def taxonomy():
			if not self.email:
				print "\nPlease enter a valid email address to let Entrez know who you are. It's *your* fault if this is not valid, and you will likely have your IP address barred from using GenBank if you don't enter one"
				Entrez.email = raw_input("")
			print "\nCreating a 'taxonomy' for your species from GenBank"
			print "\tNote: this will be saved as a file, but will not generate a constraint tree (yet...)"
			for sp in self.speciesNames:
				self.APICheck()
				self.taxonomy.append(findLineage(sp))
			print "...lineages found!"
		
		if not self.constraintMethod:
			print "It is *stronlgy* advised that you use a constraint tree with this program."
			print "\tTo supply your own constraint tree, type 'newick' and press enter."
			print "\tTo use Phyomatic to generate a constraint tree, type 'phylomatic' and press enter."
			print "\tTo generate a taxonomy for your species from GenBank, type 'taxonomy' and press enter."
			print "Otherwise, press enter to continue without a constraint tree."
			stopper = True
			while stopper:
				constraintInput = raw_input("Constraint Method: ")
				if constraintInput:
					if constraintInput == 'newic':
						stopper = newickConstraint()
					elif constraintInput == 'phylomatic':
						stopper = phylomatic()
					elif constraintInput == 'taxonomy':
						stopper = taxonomy()
					else:
						print "Sorry, I didn't get that. Please try again."
				else:
					print "...Continuing without constraint tree"
					stopper = False
		else:
			if self.constraintMethod == 'newick':
				newickConstriant()
			elif self.constraintMethod == 'phylomatic':
				phylomatic()
			elif self.constraintMethod == 'taxonomy':
				taxonomy()
			else:
				print "Error in constraint tree method name."
	
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
			for i in range(1, (len(self.genes)-1)):
				tempAlignment += self.alignment[i]
			self.alignment = tempAlignment
			for i,x in enumerate(self.speciesNames):
				self.alignment[i].id = self.speciesNames[i].replace(' ', '_')
	
	def APICheck(self):
		if self.tracker == self.spacer:
			time.sleep(self.delay)
			self.tracker = 0
		else:
			self.tracker += 1
	

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
		
		#Startup
		currentState = PhyloGenerator(args)
		
		#DNA from file
		if not len(currentState.sequences):
			print "\nDNA INPUT"
			currentState.loadDNAFile()
		
		#DNA from GenBank
		if not len(currentState.sequences):
			print "\nDNA DOWNLOAD"
			currentState.loadGenBank()
		
		#DNA Assertion
		if not len(currentState.sequences):
			print "\nNo DNA, or list of species for which to download DNA, loaded. Exiting."
			sys.exit()
		
		#DNA Checking
		"\nDNA CHECKING"
		currentState.dnaChecking()
		print "\nYou are now able to delete DNA sequences you have loaded.\nEvery time you delete a sequence, your summary statistics will be re-calculated, and displayed to you again.\n*IMPORTANT*: Sequence IDs may change once you delete a sequence."
		currentState.dnaEditing()
		
		#DNA Cleanup and renaming
		currentState.cleanUpSequences()
		currentState.renameSequences()
		
		#Alignment
		print "\nDNA ALIGNMENT"
		currentState.align()
		print "\nALIGNMENT CHECKING"
		currentState.alignmentEditing()
		currentState.concatenateSequences()
		
		#Constraint tree
		print "\nCONSTRAINT TREE"
		#VERY poorly written - fix this!!!
		currentState.getConstraint()
		
		#Phylogeny building
		currentState.phylogen()
		
		#Rate smoothing
		currentState.rateSmooth()
		
		#Output
		currentState.writeOutput()
		print "\nCongratulations! Exiting phyloGenerator."


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="phyloGenerator - phylogeny generation for ecologists.", epilog="Help at http://willpearse.github.com/phyloGenerator - written by Will Pearse")
	parser.add_argument("--version", action="store_true", help="Display version information.")
	parser.add_argument("--manual", action="store_true", help="(Attempt to) open browser and show help")
	parser.add_argument("-name", "-n", help="'Stem' name for all output files.")
	parser.add_argument("-dna", "-d", help="Unaligned DNA (in FASTA format).")
	parser.add_argument("-alignment", "-a", help="Alignment method")
	parser.add_argument("-gene", "-g", help="The genes to search for (multiple genes are comma-separated)")
	parser.add_argument("-nGenes", "-ng", help="The number of genes to search for (if fewer than suggested in 'genes' are required)")
	parser.add_argument("-species", "-s", help="Binomial names of species, each on a new line")
	parser.add_argument("-consTree", "-cT", help="Constraint tree (in newick format).")
	parser.add_argument("-consPhylomat", "-cP", help="Phylomatic's phylogeny (comma) and taxonomy")
	parser.add_argument("-email", "-e", help="Email address for GenBank searches.")
	parser.add_argument("-parameters", "-p", help="Parameter file giving detailed instructions to phyloGen.")
	parser.add_argument("-delay", help="Delay (seconds) when pausing between any internet API calls.")
	main()
