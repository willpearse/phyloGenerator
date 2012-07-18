#!/usr/bin/env python
# encoding: utf-8
"""
phyloGenerator.py
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
from scipy import percentile,mean,std
import subprocess, threading #Background process class
from Bio import AlignIO #Handle alignments
from Bio.Align import MultipleSeqAlignment #Create an alignment (for shitty MUSCLE)
import os #Remove temporary files
import re #Search for files to delete
from Bio import Phylo #Load constructed phylogeny
from Bio.Data import CodonTable #Codon positions for trimming sequences
import time #For waiting between sequence downloads
import argparse #For command line arguments
import webbrowser #Load website on request
import sys #To exit on errors
import pdb #Debugging
import warnings
maxCheck = 4
def unrootPhylomatic(tree):
	Phylo.write(tree, 'unrootingOutput.tre', 'newick')
	with open('unrootingOutput.tre') as f:
		treeText = f.read()
	if re.search('\)[a-zA-z_:0-9\.]*;', treeText):
		treeText = re.sub('\)[a-zA-z_:0-9\.]*;', ';', treeText)
		treeText = treeText[1:]
		with open('unrootingOutput.tre', 'w') as f:
			f.write(treeText)
		tree = Phylo.read('unrootingOutput.tre', 'newick')
	os.remove('unrootingOutput.tre')
	return tree

def taxonIDLookup(taxonID):
	finished = 0
	while finished <= maxCheck:
		try:
			handleDownload = Entrez.efetch(db="taxonomy", id=taxonID, retmode="xml")
			resultsDownload = Entrez.read(handleDownload)
			handleDownload.close()
			finished = maxCheck + 1
		except:
			if finished == 0:
				print "!!!Server error - retrying..."
				finished += 1
				time.sleep(10)
			if finished == maxCheck:
				print "!!!!!!Unreachable. Returning nothing."
				return(tuple())
	scientificName = resultsDownload[0]['ScientificName']
	lineage = resultsDownload[0]['Lineage'].split("; ")
	lineage.reverse()
	lineage = tuple(lineage)
	taxId = resultsDownload[0]['TaxId']
	mitoCode = resultsDownload[0]['MitoGeneticCode']['MGCName']
	return(scientificName, lineage, taxId, mitoCode)

def commonLookup(spName):
	finished = 0
	while finished <= maxCheck:
		try:
			handleSearch = Entrez.esearch(db="taxonomy", term=spName)
			resultsSearch = Entrez.read(handleSearch)
			handleSearch.close()
			finished = maxCheck + 1
		except:
			if finished == 0:
				print "!!!Server error checking", spName, " - retrying..."
				finished += 1
				time.sleep(10)
			if finished == maxCheck:
				print "!!!!!!Unreachable. Returning nothing."
				return(tuple())
	if resultsSearch['IdList']:
		return taxonIDLookup(resultsSearch['IdList'][0])
	else:
		return(tuple())

def cladeSpecies(cladeName):
	searchTerm = cladeName + '[subtree] AND species[rank]'
	finished = 0
	while finished <= maxCheck:
		try:
			handleSearch = Entrez.esearch(db="taxonomy", term=searchTerm)
			resultsSearch = Entrez.read(handleSearch)
			finished = maxCheck + 1
		except:
			if finished == 0:
				print "!!!Server error checking", cladeName, " - retrying..."
				finished += 1
				time.sleep(10)
			if finished == maxCheck:
				print "!!!!!!Unreachable. Returning nothing."
				return()
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
	finished = 0
	while finished <= maxCheck:
		try:
			handle = Entrez.esearch(db="nucleotide",term=term, usehistory=usehistory, retStart=retStart, retMax=retMax, retmode="text")
			results = Entrez.read(handle)
			handle.close()
			finished = maxCheck + 1
		except:
			if finished == 0:
				print "!!!Server error checking", term, " - retrying..."
				finished += 1
				time.sleep(1)
			if finished == maxCheck:
				print "!!!!!!Unreachable. Returning nothing."
				return()
	return results

def eFetchSeqID(seqID, rettype='gb'):
	finished = 0
	while finished <= maxCheck:
		try:
			handle = Entrez.efetch(db="nucleotide", rettype=rettype, retmode="text", id=seqID)
			results = SeqIO.read(handle,rettype)
			handle.close()
			finished = maxCheck + 1
		except:
			if finished == 0:
				print "!!!Server error - retrying..."
				finished += 1
				time.sleep(10)
			if finished == maxCheck:
				print "!!!!!!Unreachable. Returning nothing."
				return(tuple())
	return results

def eFetchESearch(eSearchOutput, rettype='gb'):
	finished = 0
	while finished <= maxCheck:
		try:
			handle = Entrez.efetch(db="nucleotide", rettype=rettype, retmode="text", webenv=eSearchOutput['WebEnv'], query_key=eSearchOutput['QueryKey'])
			results = SeqIO.read(handle, rettype)
			handle.close()
			finished = maxCheck + 1
		except:
			if finished == 0:
				print "!!!Server error - retrying..."
				finished += 1
				time.sleep(10)
			if finished == maxCheck:
				print "!!!!!!Unreachable. Returning nothing."
				return(tuple())
	return results

def eSummary(seqID):
	finished = 0
	while finished <= maxCheck:
		try:
			handle = Entrez.esummary(db="nucleotide", id=seqID, retmode="text")
			results = Entrez.read(handle)
			handle.close()
			finished = maxCheck + 1
		except:
			if finished == 0:
				print "!!!Server error - retrying..."
				finished += 1
				time.sleep(10)
			if finished == maxCheck:
				print "!!!!!!Unreachable. Returning nothing."
				return(tuple())
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
					#This next bit is necessary because you can't search for a median - a median can not be in the sample... (e.g. it could be floating point)
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
					if geneName:
						#Now we're checking to see if they haven't labelled the sequence with a gene, and instead are checking in the title for the gene
						if titleText:
							titleText += " " + geneName
						else:
							titleText = geneName
						geneName = None
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

def argsCheck(arguments, parameter, argSplit='-', paramSplit=' '):
	parameters = arguments.split(argSplit)
	for each in parameters:
		if each == parameter:
			return each.split(parameter+paramSplit)[1]
		else:
			raise RuntimeError("A match value for '" + paramter + "' was not found in the call '" + arguments + "'")

def alignSequences(seqList, sppNames, method='muscle', tempStem='temp', timeout=99999999, silent=False, nGenes=1, verbose=True):
	finalOutput = []
	output = []
	alignedSomething = False
	if method == 'everything': method = 'muscle-mafft-clustalo-prank'
	if method == 'quick': method = 'muscle-mafft-clustalo'
	for i in range(nGenes):
		if verbose: print "...aligning gene no.", i+1
		geneOutput = []
		seqs = [x[i] for x in seqList]
		if 'muscle' in method:
			print "......with MUSCLE"
			inputFile = tempStem + '.fasta'
			outputFile = tempStem + 'Out.fasta'
			commandLine = 'muscle -in ' + inputFile + " -out " + outputFile
			SeqIO.write(seqs, inputFile, "fasta")
			pipe = TerminationPipe(commandLine, timeout)
			pipe.run(silent=silent)
			os.remove(inputFile)
			if not pipe.failure:
				#Need to put the sequences back in order (...)
				try:
					newSeqs = AlignIO.read(outputFile, 'fasta')
				except:
					raise RuntimeError("MUSCLE unable to run: check your input sequences please!")
				sortedAligns = []
				for oldSeq in seqs:
					for alignSeq in newSeqs:
						if alignSeq.name == oldSeq.name or alignSeq.name == oldSeq.id:
							sortedAligns.append(alignSeq)
				geneOutput.append(MultipleSeqAlignment(sortedAligns))
				os.remove(outputFile)
				alignedSomething = True
			else:
				raise RuntimeError("MUSCLE alignment not complete in time allowed")
		
		if 'mafft' in method:
			print "......with MAFFT"
			inputFile = tempStem + '.fasta'
			outputFile = tempStem + 'Out.fasta'
			commandLine = 'mafft --auto ' + inputFile + " > " + outputFile
			SeqIO.write(seqs, inputFile, "fasta")
			pipe = TerminationPipe(commandLine, timeout)
			pipe.run(silent=silent)
			os.remove(inputFile)
			if not pipe.failure:
				try:
					geneOutput.append(AlignIO.read(outputFile, 'fasta'))
				except:
					raise RuntimeError("MAFFT unable to run: check your input sequences please!")
				os.remove(outputFile)
				alignedSomething = True
			else:
				raise RuntimeError("Mafft alignment not complete in time allowed")
		
		if 'clustalo' in method:
			print "......with Clustal-Omega"
			inputFile = tempStem + '.fasta'
			outputFile = tempStem + 'Out.fasta'
			commandLine = 'clustalo -i ' + inputFile + " -o " + outputFile + " -v"
			#Remove missing sequences
			outSeqs = []
			for seq in seqs:
				if seq.seq.tostring() != '':
					outSeqs.append(seq)
			SeqIO.write(outSeqs, inputFile, "fasta")
			pipe = TerminationPipe(commandLine, timeout)
			pipe.run(silent=silent)
			os.remove(inputFile)
			if not pipe.failure:
				try:
					clustalAlign = AlignIO.read(outputFile, 'fasta')
					alignNames = [x.id for x in clustalAlign]
					finalList = []
					x = 0
					fillerSeq = SeqRecord(Seq('-' * clustalAlign.get_alignment_length()))
					for i,seq in enumerate(seqs):
						if seq.id in alignNames:
							finalList.append(clustalAlign[x])
							x += 1
						else:
							finalList.append(fillerSeq)
							finalList[-1].id = sppNames[i] + "MadeUp"
							finalList[-1].name = sppNames[i] + "MadeUp"
							finalList[-1].description = "Made up by phyloGenerator"
					os.remove(outputFile)
					
					geneOutput.append(MultipleSeqAlignment(finalList))
					alignedSomething = True
				except:
					raise RuntimeError("Clustal-Omega unable to run: check your input sequences please!")
			else:
				raise RuntimeError("Clustal-o alignment not complete in time allowed")
		
		if 'prank' in method:
			print "......with Prank"
			inputFile = tempStem + '.fasta'
			outputFile = tempStem + 'Out.fasta'
			commandLine = 'prank -d=' + inputFile + " -o=" + outputFile
			SeqIO.write(seqs, inputFile, "fasta")
			pipe = TerminationPipe(commandLine, timeout)
			pipe.run(silent=silent)
			os.remove(inputFile)
			if not pipe.failure:
				try:
					geneOutput.append(AlignIO.read(outputFile+".2.fas", 'fasta'))
				except:
					raise RuntimeError("PRANK unable to run: check your input sequences please!")
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
			output['quantileLengths'].append(percentile(each, [5, 25, 50, 75, 95]))
			output['maxLength'].append(max(each))
			output['minLength'].append(min(each))
			output['upperQuantile'].append([False] * len(each))
			output['lowerQuantile'].append([False] * len(each))
			for k in range(len(each)):
				if each[k] <= output['quantileLengths'][i][0]:
					output['lowerQuantile'][i][k] = True
				elif each[k] >= output['quantileLengths'][i][4]:
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
	if maxInputName < len("Input Name"): maxInputName = len("Input Name")
	maxInputName += 1
	#Setup geneNames lengths
	geneNameLengths = [len(x) for x in geneNames]
	for i,each in enumerate(geneNameLengths):
		if each < 6:
			geneNameLengths[i] = 6
	
	#Print out details, doing things differently if we haven't analysed them
	if(seqDetect):
		#if seqDetect['tolerable']:
		#	print "Sequence lengths within specified tolerance"
		#else:
		#	print "Sequence lengths *outside* specified tolerance"
		print "\nSequence summary:"
		print "'^^^' and '___' denote particularly long or short sequences\n"
		header = "Sp. ID " + "Input name".ljust(maxInputName)
		for i, each in enumerate(geneNames):
			header += each.ljust(geneNameLengths[i]) + "     "
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
						stars.append("   ")
				else:
					#What *does* happen to empty sequences?
					stars.append("   ")
			
			row = str(i).ljust(len("Sp. ID ")) + str(speciesNames[i]).ljust(maxInputName)
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
	print ""

def alignmentDisplay(alignments, alignMethods, geneNames, alignDetect=None):
	assert len(alignments)==len(geneNames)
	assert len(alignments[0])==len(alignMethods)
	for alignNo, alignList in enumerate(alignments):
		print "\nGene:", geneNames[alignNo]
		if alignDetect:
			print "ID", "Alignment".ljust(12), "Length".ljust(10), "Med. Gaps".ljust(15), "SD Gaps".ljust(15), "Min-Max Gaps".ljust(15), "Med. Gap Frac.".ljust(15), "Max Gap Frac.".ljust(15)
			for i in range(len(alignList)):
				print  '{ID:3}{alignment:<13}{length:<11}{medGaps:<16.1f}{sdGaps:<16.2f}{minMaxGaps:16}{medGapsFrac:<16.2f}{minMaxGapsFrac:15}'.format(ID=str(i), alignment=alignMethods[i], length=alignDetect['length'][alignNo][i], medGaps=alignDetect['noGaps'][alignNo]['median'][i], sdGaps=alignDetect['noGaps'][alignNo]['sd'][i], minMaxGaps=str(str(round(alignDetect['noGaps'][alignNo]['min'][i],3))+" - "+str(round(alignDetect['noGaps'][alignNo]['max'][i],3))), medGapsFrac=alignDetect['gapFraction'][alignNo]['median'][i], minMaxGapsFrac=str(str(round(alignDetect['gapFraction'][alignNo]['min'][i],3))+"-"+str(round(alignDetect['gapFraction'][alignNo]['max'][i],3))))
		else:
			print "ID", "Alignment", "Length"
			for i in range(len(alignList)):
				print str(i).ljust(len("ID")), alignMethods[i].ljust(len("Alignment")), str(alignList[i].get_alignment_length()).ljust(len("Length"))

def RAxML(alignment, method='localVersion', tempStem='temp', outgroup=None, timeout=999999999, cladeList=None,	DNAmodel='GTR+G', partitions=None, constraint=None, cleanup=True, runNow=True, startingTree=None):
	inputFile = tempStem + 'In.phylip'
	outputFile = tempStem + 'Out'
	fileLine = ' -s ' + inputFile + ' -n ' + outputFile
	options = ' -p ' + str(random.randint(0,10000000))
	#Outgroup(s), assuming they're in the right format for RAxML!
	if outgroup:
		options += ' -o ' + outgroup
	AlignIO.write(alignment, inputFile, "phylip-relaxed")
	if 'startingOnly' in method:
		algorithm = ''
		DNAmodel = ''
		options += ' -y'
	else:
		#Restarts
		if 'restarts' in method:
			temp = method.split('-')
			for each in temp:
				if 'restarts' in each:
					number = each.replace('restarts=', '')
					options += ' -N ' + number
					break
		#Algorithm
		if 'integratedBootstrap' in method:
			if 'restarts' in method:
				raise RuntimeError("RAxML cannot do the accelerated bootstrap multiple times (at least not the way I'm using it).")
			temp = method.split('-')
			for each in temp:
				if 'integratedBootstrap' in each:
					number = each.replace('integratedBootstrap=', '')
					options += ' -N ' + number
					break
			options += ' -x ' + str(random.randint(0,10000000))
			algorithm = ' -f a'
		else:
			algorithm = ' -f d'
		#Starting tree
		if startingTree:
			Phylo.write(startingTree, tempStem + 'startingTree.tre', 'newick')
			options += ' -t ' + tempStem + 'startingTree.tre'
	#RAxML Flavour
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
	#Partitions
	if partitions:
		with open(tempStem+"_partitions.txt", 'w') as f:
			for i in range(0, len(partitions)-1):
				f.write("DNA, position" + str(partitions[i]) + " = " + str(partitions[i]+1) + "-" + str(partitions[i+1]) + "\n")
		options += " -q " + tempStem + "_partitions.txt"
	
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
		if 'startingOnly' in method:
			tree = Phylo.read('RAxML_parsimonyTree.' + outputFile, "newick")
		elif 'integratedBootstrap' in method:
			tree = Phylo.read('RAxML_bipartitions.' + outputFile, "newick")
		else:
			tree = Phylo.read('RAxML_bestTree.' + outputFile, "newick")
		if cleanup:
			if startingTree:
				os.remove(tempStem + 'startingTree.tre')
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

def BEAST(alignment, method='GTR+GAMMA', tempStem='temp', timeout=999999999, constraint=None, cleanup=False, runNow=True, chainLength=1000000, logRate=1000, screenRate=1000, overwrite=True, burnin=0.1, restart=None):
	completeConstraint = False
	agedConstraint = False
	with open(tempStem+"_BEAST.xml", 'w') as f:
		f.write('<?xml version="1.0" standalone="yes"?>\n')
		f.write('<beast>\n')
		f.write('	<!-- The list of taxa analyse (can also include dates/ages).				 -->\n')
		f.write('	<taxa id="taxa">\n')
		if type(alignment) is list:
			for each in alignment[0]:
				f.write('		<taxon id="'+each.id+'"/>\n')
		else:
			for each in alignment:
				f.write('		<taxon id="'+each.id+'"/>\n')
		f.write('	</taxa>\n')
		if constraint:
			clades = []
			cladesNames = []
			cladesAges = []
			agedNames = []
			ages = []
			for x,clade in enumerate(constraint.find_clades()):
				temp = [x.name for x in clade.get_terminals()]
				if len(temp) > 1:
					clades.append(temp)
					cladesNames.append(clade.name)
					cladesAges.append(clade.branch_length)
			if len(set([item for sublist in clades for item in sublist])) == len(constraint.get_terminals()):
				completeConstraint = True
			for i,clade in enumerate(clades):
				if cladesNames[i]:
					ages.append(cladesAges[i])
					agedNames.append("cladeNo"+str(i))
					agedConstraint = True
				f.write('	<taxa id="cladeNo' + str(i) + '">\n')
				for sp in clade:
					f.write('		<taxon idref="' + sp + '"/>\n')
				f.write('	</taxa>')
		f.write('	<!-- The sequence alignment (each sequence refers to a taxon above).		 -->\n')
		if type(alignment) is list:
			for i,align in enumerate(alignment):
				f.write('	<alignment id="alignment'+str(i)+'" dataType="nucleotide">\n')
				for seq in align:
					f.write('		<sequence>\n')
					f.write('			<taxon idref="'+seq.id+'"/>\n')
					f.write(seq.seq.tostring()+'\n')
					f.write('		</sequence>\n')
				f.write('	</alignment>\n')
		else:
			f.write('	<alignment id="alignment" dataType="nucleotide">\n')
			for each in alignment:
				f.write('		<sequence>\n')
				f.write('			<taxon idref="'+each.id+'"/>\n')
				f.write(each.seq.tostring()+'\n')
				f.write('		</sequence>\n')
			f.write('	</alignment>\n')
		f.write('	<!-- The unique patterns from 1 to end										 -->\n')
		if type(alignment) is list:
			for i,align in enumerate(alignment):
				f.write('	<patterns id="patterns' + str(i) + '" from="1">\n')
				f.write('		<alignment idref="alignment' + str(i) + '"/>\n')
				f.write('	</patterns>\n')
		else:
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
			Phylo.write(constraint, 'tempStem'+'_TREE.tre', 'newick')
			with open('tempStem'+'_TREE.tre') as tF:
				treeFormat = tF.readlines()[0]
				f.write('	<newick id="startingTree">')
				f.write('		'+treeFormat)
				f.write('	</newick>')
			os.remove('tempStem'+'_TREE.tre')
		else:
			f.write('	<!-- Generate a random starting tree under the coalescent process			 -->\n')
			f.write('	<coalescentTree id="startingTree" rootHeight="0.092">\n')
			if constraint:
				for i,clade in enumerate(clades):
					f.write('	<constrainedTaxa>\n')
					f.write('		<tmrca monophylletic="true">\n')
					f.write('			<taxa idref="cladeNo' + str(i) + '">\n')
					f.write('		</tmrca>\n')
			f.write('		<taxa idref="taxa"/>\n')
			if constraint:
				f.write('	</constrainedTaxa>\n')
			f.write('		<constantSize idref="initialDemo"/>\n')
			f.write('	</coalescentTree>\n')
		f.write('	<!-- Generate a tree model													 -->\n')
		f.write('	<treeModel id="treeModel">')
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
		if constraint:
			for i,clade in enumerate(clades):
				f.write('	<!-- Taxon Sets                                                              -->\n')
				f.write('	<tmrcaStatistic id="tmrca(cladeNo' + str(i) + ')" includeStem="false">\n')
				f.write('			<mrca>\n')
				f.write('				<taxa idref="cladeNo' + str(i) + '"/>\n')
				f.write('			</mrca>\n')
				f.write('		<treeModel idref="treeModel"/>\n')
				f.write('	</tmrcaStatistic>\n')
				f.write('	<monophylyStatistic id="monophyly(' + str(i) + ')">\n')
				f.write('		<mrca>\n')
				f.write('			<taxa idref="cladeNo' + str(i) + '"/>\n')
				f.write('		</mrca>\n')
				f.write('<treeModel idref="treeModel"/>\n')
				f.write('</monophylyStatistic>\n')
		f.write('	<!-- Generate a speciation likelihood for Yule or Birth Death				 -->\n')
		f.write('	<speciationLikelihood id="speciation">\n')
		f.write('		<model>\n')
		f.write('			<yuleModel idref="yule"/>\n')
		f.write('		</model>\n')
		f.write('		<speciesTree>\n')
		f.write('			<treeModel idref="treeModel"/>\n')
		f.write('		</speciesTree>\n')
		f.write('	</speciationLikelihood>\n')
		f.write('	<!-- Uncorrelated relaxed clock models (Drummond, Ho, Phillips & Rambaut, 2006) -->\n')
		if type(alignment) is list:
			for i,align in enumerate(alignment):
				f.write('	<discretizedBranchRates id="branchRates' + str(i) + '">\n')
				f.write('		<treeModel idref="treeModel"/>\n')
				f.write('		<distribution>\n')
				f.write('			<logNormalDistributionModel meanInRealSpace="true">\n')
				f.write('				<mean>\n')
				f.write('					<parameter id="ucld.mean' + str(i) + '" value="1.0"/>\n')
				f.write('				</mean>\n')
				f.write('				<stdev>\n')
				f.write('					<parameter id="ucld.stdev' + str(i) + '" value="0.3333333333333333" lower="0.0" upper="Infinity"/>\n')
				f.write('				</stdev>\n')
				f.write('			</logNormalDistributionModel>\n')
				f.write('		</distribution>\n')
				f.write('		<rateCategories>\n')
				f.write('			<parameter id="branchRates' + str(i) + '.categories" dimension="18"/>\n')
				f.write('		</rateCategories>\n')
				f.write('	</discretizedBranchRates>\n')
				
				f.write('	<rateStatistic id="meanRate' + str(i) + '" name="meanRate" mode="mean" internal="true" external="true">\n')
				f.write('		<treeModel idref="treeModel"/>\n')
				f.write('		<discretizedBranchRates idref="branchRates' + str(i) + '"/>\n')
				f.write('	</rateStatistic>\n')
				
				f.write('	<rateStatistic id="coefficientOfVariation' + str(i) + '" name="coefficientOfVariation" mode="coefficientOfVariation" internal="true" external="true">\n')
				f.write('		<treeModel idref="treeModel"/>\n')
				f.write('		<discretizedBranchRates idref="branchRates' + str(i) + '"/>\n')
				f.write('	</rateStatistic>\n')
				
				f.write('	<rateCovarianceStatistic id="covariance' + str(i) + '" name="covariance">\n')
				f.write('		<treeModel idref="treeModel"/>\n')
				f.write('		<discretizedBranchRates idref="branchRates' + str(i) + '"/>\n')
				f.write('	</rateCovarianceStatistic>\n')
		else:
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
		if type(alignment) is list:
			for i,align in enumerate(alignment):
				if 'GTR' in method:
					f.write('	<!-- The general time reversible (GTR) substitution model					 -->\n')
					f.write('	<gtrModel id="gtr' + str(i) + '">\n')
					f.write('		<frequencies>\n')
					f.write('			<frequencyModel dataType="nucleotide">\n')
					f.write('				<frequencies>\n')
					if 'base=empirical' in method:
						f.write('					<parameter id="frequencies' + str(i) + '" dimension="4"/>\n')
					else:
						f.write('					<parameter id="frequencies' + str(i) + '" value="0.25 0.25 0.25 0.25"/>\n')
					f.write('				</frequencies>\n')
					f.write('			</frequencyModel>\n')
					f.write('		</frequencies>\n')
					f.write('		<rateAC>\n')
					f.write('			<parameter id="ac' + str(i) + '" value="1.0" lower="0.0" upper="Infinity"/>\n')
					f.write('		</rateAC>\n')
					f.write('		<rateAG>\n')
					f.write('			<parameter id="ag' + str(i) + '" value="1.0" lower="0.0" upper="Infinity"/>\n')
					f.write('		</rateAG>\n')
					f.write('		<rateAT>\n')
					f.write('			<parameter id="at' + str(i) + '" value="1.0" lower="0.0" upper="Infinity"/>\n')
					f.write('		</rateAT>\n')
					f.write('		<rateCG>\n')
					f.write('			<parameter id="cg' + str(i) + '" value="1.0" lower="0.0" upper="Infinity"/>\n')
					f.write('		</rateCG>\n')
					f.write('		<rateGT>\n')
					f.write('			<parameter id="gt' + str(i) + '" value="1.0" lower="0.0" upper="Infinity"/>\n')
					f.write('		</rateGT>\n')
					f.write('	</gtrModel>\n')
				elif 'HKY' in method:
					f.write('	<!-- The HKY substitution model (Hasegawa, Kishino & Yano, 1985)			 -->')
					f.write('	<HKYModel id="hky' + str(i) + '">')
					f.write('		<frequencies>')
					f.write('			<frequencyModel dataType="nucleotide">')
					f.write('				<frequencies>')
					f.write('					<parameter id="frequencies' + str(i) + '" value="0.25 0.25 0.25 0.25"/>')
					f.write('				</frequencies>')
					f.write('			</frequencyModel>')
					f.write('		</frequencies>')
					f.write('		<kappa>')
					f.write('			<parameter id="kappa' + str(i) + '" value="2.0" lower="0.0" upper="Infinity"/>')
					f.write('		</kappa>')
					f.write('	</HKYModel>')
				else:
					raise RuntimeError("No valid DNA substituion model specified for BEAST.")
		else:
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
				f.write('	<!-- The HKY substitution model (Hasegawa, Kishino & Yano, 1985)			 -->')
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
		if type(alignment) is list:
			for i,align in enumerate(alignment):
				f.write('	<!-- site model																 -->\n')
				f.write('	<siteModel id="siteModel' + str(i) + '">\n')
				f.write('		<substitutionModel>\n')
				if 'GTR' in method:
					f.write('			<gtrModel idref="gtr' + str(i) + '"/>\n')
					if 'GAMMA' in method:
						f.write('		<gammaShape gammaCategories="4">')
						f.write('			<parameter id="alpha' + str(i) + '" value="0.5" lower="0.0" upper="1000.0"/>')
						f.write('		</gammaShape>')
				elif 'HKY' in method:
					f.write('			<HKYModel idref="hky' + str(i) + '"/>')
					if 'GAMMA' in method:
						f.write('		<gammaShape gammaCategories="4">')
						f.write('			<parameter id="alpha' + str(i) + '" value="0.5" lower="0.0" upper="1000.0"/>')
						f.write('		</gammaShape>')
				else:
					raise RuntimeError("No valid DNA substituion model specified for BEAST.")
				f.write('		</substitutionModel>\n')
			f.write('	</siteModel>\n')
		else:
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
		if type(alignment) is list:
			for i,align in enumerate(alignment):
				f.write('	<treeLikelihood id="treeLikelihood' + str(i) + '" useAmbiguities="false">\n')
				f.write('		<patterns idref="patterns' + str(i) + '"/>\n')
				f.write('		<treeModel idref="treeModel"/>\n')
				f.write('		<siteModel idref="siteModel' + str(i) + '"/>\n')
				f.write('		<discretizedBranchRates idref="branchRates' + str(i) + '"/>\n')
				f.write('	</treeLikelihood>\n')
		else:
			f.write('	<treeLikelihood id="treeLikelihood" useAmbiguities="false">\n')
			f.write('		<patterns idref="patterns"/>\n')
			f.write('		<treeModel idref="treeModel"/>\n')
			f.write('		<siteModel idref="siteModel"/>\n')
			f.write('		<discretizedBranchRates idref="branchRates"/>\n')
			f.write('	</treeLikelihood>\n')
		f.write('	<!-- Define operators														 -->\n')
		f.write('	<operators id="operators">\n')
		if type(alignment) is list:
			for i,align in enumerate(alignment):
				if 'GTR' in method:
					f.write('		<scaleOperator scaleFactor="0.75" weight="3">\n')
					f.write('			<parameter idref="ucld.stdev' + str(i) + '"/>\n')
					f.write('		</scaleOperator>\n')
					f.write('		<scaleOperator scaleFactor="0.75" weight="0.1">\n')
					f.write('			<parameter idref="ac' + str(i) + '"/>\n')
					f.write('		</scaleOperator>\n')
					f.write('		<scaleOperator scaleFactor="0.75" weight="0.1">\n')
					f.write('			<parameter idref="ag' + str(i) + '"/>\n')
					f.write('		</scaleOperator>\n')
					f.write('		<scaleOperator scaleFactor="0.75" weight="0.1">\n')
					f.write('			<parameter idref="at' + str(i) + '"/>\n')
					f.write('		</scaleOperator>\n')
					f.write('		<scaleOperator scaleFactor="0.75" weight="0.1">\n')
					f.write('			<parameter idref="cg' + str(i) + '"/>\n')
					f.write('		</scaleOperator>\n')
					f.write('		<scaleOperator scaleFactor="0.75" weight="0.1">\n')
					f.write('			<parameter idref="gt' + str(i) + '"/>\n')
					f.write('		</scaleOperator>\n')
				elif 'HKY' in method:
					f.write('		<deltaExchange delta="0.01" weight="0.1">\n')
					f.write('			<parameter idref="frequencies' + str(i) + '"/>\n')
					f.write('		</deltaExchange>\n')
				else:
					raise RuntimeError("No valid DNA substituion model specified for BEAST.")
		else:
			f.write('		<scaleOperator scaleFactor="0.75" weight="3">\n')
			f.write('			<parameter idref="ucld.stdev"/>\n')
			f.write('		</scaleOperator>\n')
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
		if type(alignment) is list:
			for i,aling in enumerate(alignment):
				f.write('		<swapOperator size="1" weight="10" autoOptimize="false">\n')
				f.write('			<parameter idref="branchRates' + str(i) + '.categories"/>\n')
				f.write('		</swapOperator>\n')
				f.write('		<randomWalkIntegerOperator windowSize="1" weight="10">\n')
				f.write('			<parameter idref="branchRates' + str(i) + '.categories"/>\n')
				f.write('		</randomWalkIntegerOperator>\n')
				f.write('		<uniformIntegerOperator weight="10">\n')
				f.write('			<parameter idref="branchRates' + str(i) + '.categories"/>\n')
				f.write('		</uniformIntegerOperator>\n')
		else:
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
		if agedConstraint:
			for age,clade in zip(ages, agedNames):
				f.write('				<normalPrior mean="'+ str(age) + '" stdev="1.0">\n')
				f.write('					<statistic idref="tmrca(' + str(clade) + ')"/>\n')
				f.write('				</normalPrior>\n')
		if 'GTR' in method:
			if type (alignment) is list:
				for i,align in enumerate(alignment):
					f.write('				<gammaPrior shape="0.05" scale="10.0" offset="0.0">\n')
					f.write('					<parameter idref="ac' + str(i) + '"/>\n')
					f.write('				</gammaPrior>\n')
					f.write('				<gammaPrior shape="0.05" scale="20.0" offset="0.0">\n')
					f.write('					<parameter idref="ag' + str(i) + '"/>\n')
					f.write('				</gammaPrior>\n')
					f.write('				<gammaPrior shape="0.05" scale="10.0" offset="0.0">\n')
					f.write('					<parameter idref="at' + str(i) + '"/>\n')
					f.write('				</gammaPrior>\n')
					f.write('				<gammaPrior shape="0.05" scale="10.0" offset="0.0">\n')
					f.write('					<parameter idref="cg' + str(i) + '"/>\n')
					f.write('				</gammaPrior>\n')
					f.write('				<gammaPrior shape="0.05" scale="10.0" offset="0.0">\n')
					f.write('					<parameter idref="gt' + str(i) + '"/>\n')
					f.write('				</gammaPrior>\n')
					f.write('				<exponentialPrior mean="0.3333333333333333" offset="0.0">\n')
					f.write('					<parameter idref="ucld.stdev' + str(i) + '"/>\n')
					f.write('				</exponentialPrior>\n')
			else:
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
		if type(alignment) is list:
			for i,align in enumerate(alignment):
				f.write('				<treeLikelihood idref="treeLikelihood' + str(i) + '"/>\n')
		else:
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
		if type(alignment) is list:
				for i,align in enumerate(alignment):
					f.write('			<column label="ucld.mean' + str(i) + '" sf="6" width="12">\n')
					f.write('				<parameter idref="ucld.mean' + str(i) + '"/>\n')
					f.write('			</column>\n')
		else:
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
			if type(alignment) is list:
				for i,align in enumerate(alignment):
					f.write('			<parameter idref="ac' + str(i) + '"/>\n')
					f.write('			<parameter idref="ag' + str(i) + '"/>\n')
					f.write('			<parameter idref="at' + str(i) + '"/>\n')
					f.write('			<parameter idref="cg' + str(i) + '"/>\n')
					f.write('			<parameter idref="gt' + str(i) + '"/>\n')
					f.write('			<parameter idref="frequencies' + str(i) + '"/>\n')
			else:
				f.write('			<parameter idref="ac"/>\n')
				f.write('			<parameter idref="ag"/>\n')
				f.write('			<parameter idref="at"/>\n')
				f.write('			<parameter idref="cg"/>\n')
				f.write('			<parameter idref="gt"/>\n')
				f.write('			<parameter idref="frequencies"/>\n')
		elif 'HKY' in method:
			if type(alignment) is list:
				for i,align in enumerate(alignment):
					f.write('			<parameter idref="kappa' + str(i) + '"/>')
			else:
				f.write('			<parameter idref="kappa"/>')
		else:
			raise RuntimeError("No valid DNA substituion model specified for BEAST.")
		if 'GAMMA' in method:
			f.write('			<parameter idref="alpha"/>')
		if type(alignment) is list:
			for i,align in enumerate(alignment):
				f.write('			<parameter idref="ucld.mean' + str(i) + '"/>\n')
				f.write('			<parameter idref="ucld.stdev' + str(i) + '"/>\n')
				f.write('			<rateStatistic idref="meanRate' + str(i) + '"/>\n')
				f.write('			<rateStatistic idref="coefficientOfVariation' + str(i) + '"/>\n')
				f.write('			<rateCovarianceStatistic idref="covariance' + str(i) + '"/>\n')
		else:
			f.write('			<parameter idref="ucld.mean"/>\n')
			f.write('			<parameter idref="ucld.stdev"/>\n')
			f.write('			<rateStatistic idref="meanRate"/>\n')
			f.write('			<rateStatistic idref="coefficientOfVariation"/>\n')
			f.write('			<rateCovarianceStatistic idref="covariance"/>\n')
		if type(alignment) is list:
				for i,align in enumerate(alignment):
					f.write('			<treeLikelihood idref="treeLikelihood' + str(i) + '"/>\n')
		else:
			f.write('			<treeLikelihood idref="treeLikelihood"/>\n')
		f.write('			<speciationLikelihood idref="speciation"/>\n')
		f.write('		</log>\n')
		f.write('		<!-- write tree log to file													 -->\n')
		f.write('		<logTree id="treeFileLog" logEvery="' + str(logRate) + '" nexusFormat="true" fileName="'+tempStem+'.trees" sortTranslationTable="true">\n')
		f.write('			<treeModel idref="treeModel"/>\n')
		if type(alignment) is list:
			for i,align in enumerate(alignment):
				f.write('			<discretizedBranchRates idref="branchRates' + str(i) + '"/>\n')
		else:
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
	
	if overwrite:
		commandLine = "beast -overwrite " + tempStem + "_BEAST.xml"
	else:
		commandLine = "beast -overwrite " + tempStem + "_BEAST.xml"
	
	if not runNow:
		return commandLine
	
	pipe = TerminationPipe(commandLine, timeout)
	if sys.platform == "win32":
                pipe.run(silent=False, changeDir=False)
        else:
                pipe.run(silent=False)
	if not pipe.failure:
		print "...removing burn-in of ", str(burnin*100), "%..."
		burnin = int(burnin * (chainLength / logRate))
		commandLine = 'treeannotator -burnin ' + str(burnin) + ' -heights median ' + tempStem + '.trees ' + tempStem + 'Final.tre'
		pipeAnotate = TerminationPipe(commandLine, timeout)
		pipeAnotate.run(changeDir=False)
		if not pipeAnotate.failure:
			if cleanup:
				os.remove(tempStem + "_BEAST.xml")
				os.remove(tempStem + ".trees")
				os.remove(tempStem + ".log")
				os.remove("mcmc.operators")
				return tempStem + "Final.tre"
			else:
				os.remove("mcmc.operators")
				return tempStem + "Final.tre", tempStem + "_BEAST.xml", tempStem + ".trees", tempStem + ".log"
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

def findGeneInSeq(seq, gene, trimSeq=False, DNAtype='Standard', gapType='-', verbose=False):
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
			if verbose: warnings.warn("Gene '" + gene + "' not found in sequence")
			return ()
	else:
		raise RuntimeError('No sequence features found: are you using a GenBank record?')

def cleanSequenceWrapper(seq, gene, DNAtype='Standard', gapType='-'):
	output = findGeneInSeq(seq, gene, DNAtype=DNAtype, gapType=gapType, verbose=False)
	if len(output) == len(seq):
		output = trimSequence(seq, DNAtype=DNAtype, gapType=gapType)
	return output

def rateSmooth(phylo, method='PATHd8', nodes=tuple(), sequenceLength=int(), tempStem='temp', timeout=999999):
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
				species = map(lambda x: x.name, phylo.get_terminals())
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
		output = []
		for i,gene in enumerate(align):
			geneOutput = []
			for j,method in enumerate(gene):
				AlignIO.write(method, tempStem+"Input.fasta", "fasta")
				fileLine = " -in " + tempStem + "Input.fasta -out " + tempStem + "Output.fasta -fasta"
				trimalVersion = "trimal"
				commandLine = trimalVersion + fileLine + options
				if timeout:
					pipe = TerminationPipe(commandLine, timeout)
					pipe.run()
					if not pipe.failure:
						geneOutput.append(AlignIO.read(tempStem + "Output.fasta", "fasta"))
						os.remove(tempStem + "Output.fasta")
						os.remove(tempStem + "Input.fasta")
					else:
						raise RuntimeError("Either trimAl failed, or it ran out of time")
			output.append(geneOutput)
		return output
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

def metal(alignList, tempStem='tempMetal', timeout=100):
	#Write out alignments
	alignLocations = []
	for i,gene in enumerate(alignList):
		for j,method in enumerate(gene):
			AlignIO.write(method, str(i)+"_"+str(j)+"_"+tempStem+".fasta", "fasta")
			alignLocations.append(str(i)+"_"+str(j)+"_"+tempStem+".fasta")
	
	#Use each alignment and print it out
	distances = []
	for i,currAlign in enumerate(alignLocations):
		currDist = []
		for j,nextAlign in enumerate(alignLocations[(i+1):]):
			pipe = TerminationPipe("metal --ignore-names "+currAlign+" "+nextAlign, timeout=timeout)
			pipe.run()
			temp = re.search('[0-9]*\ /\ [0-9]*', pipe.output[0]).group()
			temp = 'float(' + temp + ')'
			currDist.append(eval(temp))
		if len(currDist) > 0: distances.append(currDist)
	for i,location in enumerate(alignLocations):
		os.remove(location)
	return distances

class TerminationPipe(object):
	#Background process class
	def __init__(self, cmd, timeout, silent=True):
		self.cmd = cmd
		self.timeout = timeout
		self.process = None
		self.output = None
		self.failure = False
		self.stderr = 'EMPTY'
		self.stdout = 'EMPTY'
		self.silent = silent
	
	def run(self, silent=None, changeDir=True):
		def silentTarget():
			if sys.platform == 'win32':
                                if changeDir:
                                        self.process = subprocess.Popen("requires\\" + self.cmd, stdout=subprocess.PIPE, shell=True, stderr=subprocess.PIPE)
                                else:
                                        self.process = subprocess.Popen(self.cmd, stdout=subprocess.PIPE, shell=True, stderr=subprocess.PIPE)
			else:
                                if changeDir:
                                        self.process = subprocess.Popen("./requires/" + self.cmd, stdout=subprocess.PIPE,shell=True, stderr=subprocess.PIPE)
                                else:
                                        self.process = subprocess.Popen(self.cmd, stdout=subprocess.PIPE,shell=True, stderr=subprocess.PIPE)
			self.output = self.process.communicate()
		
		def loudTarget():
			if sys.platform == 'win32':
                                if changeDir:
                                        self.process = subprocess.Popen("requires\\" + self.cmd, shell=False)
                                else:
                                        self.process = subprocess.Popen(self.cmd, shell=False)
			else:
                                if changeDir:
                                        self.process = subprocess.Popen("./requires/" + self.cmd, shell=False)
                                else:
                                        self.process = subprocess.Popen(self.cmd, shell=False)
			self.output=self.process.communicate()
		
		if silent: self.silent = silent
		if self.silent:
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
		self.phylogenyMerged = False
		self.alignment = []
		self.alignmentUnpartitioned = []
		self.smoothPhylogeny = []
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
		self.phylogenyMethods = ''
		self.initialSeqChoice = 'medianLength'
		self.replaceSeqChoice = 'targetLength'
		self.uniqueTaxonomy = []
		self.tracker = 0
		self.rateSmoothMethods = ''
		self.allAlignments = False
		self.partitions = None
		self.alignRF = []
		self.mergedSpp = []
		self.constraintFile = False
		self.smoothPhylogeny = False
		self.smoothPhylogenyMerged = False
		#Backup working directory
		self.oldDirectory = os.getcwd()
		
		#Stem name
		if args.name:
			self.stem = args.name
			print "\nUsing stem name '" + args.name + "'..."
		else:
			print "\nPlease input a 'stem' name to act as a prefix to all output (e.g., 'stemName_phylogeny.tre')\n"
			self.stem = raw_input("Stem name: ")
		
		#Working directory
		if args.wd:
			self.workingDirectory = args.wd
			print "\nUsing working directory '" + args.wd + "' ..."
		else:
			print "\nPlease input a working directory for all your output\n"
			self.workingDirectory = raw_input("Working directory: ")
		
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
			print "\nUsing gene(s)", self.genes
			self.codonModels = ['Standard'] * len(self.genes)
		else:
			print "Please enter the gene(s) you want to use (e.g., 'COI' for cytochrome oxidase one')"
			print "Each gene on a separate line, and an empty line to continue\n"
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
		
		if args.phylogen:
			self.phylogenyMethods = args.phylogen
		
		#Constraint tree
		self.constraintMethod = ''
		#Phylomatic
		if args.consPhylomat:
			self.phylomaticPhylogeny, self.phylomaticTaxonomy = args.consPhylomat.split(",")
			self.phylomaticPhylogeny = self.phylomaticPhylogeny
			self.phylomaticTaxonomy = self.phylomaticTaxonomy
			self.constraintMethod = 'phylomatic'
		else:
			self.phylomaticPhylogeny = None
			self.phylomaticTaxonomy = None
		
		self.constraintRFs = []
		
		#Pre-supplied
		if args.consTree:
			self.constraintFile = args.consTree
			self.constraintMethod = 'newick'
			self.constraint = Phylo.read(self.constraintFile, 'newick')
		
		#Options file
		if args.options:
			options = []
			try:
				with open(args.options) as f:
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
				print "...DNA loaded"
				if not self.genes:
					print "Please enter the name of the gene you're using below\n"
					self.genes.append(raw_input("Gene name: "))
					self.nGenes = 1
				self.codonModels.append('Standard')
			except IOError:
				print "\nDNA sequence file not found. Exiting..."
				sys.exit()
		else:
			locker = True
			print "\nIf you already have DNA sequences in a FASTA file, please enter its location"
			print "Otherwise, hit enter to continue\n"
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
				sys.exit()
		else:
			locker = True
			aborted = False
			print "\nPlease enter the location of the list of species for which you want to build a phylogeny"
			print "Each species must be on a new line\n"
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
			print "\n", len(self.speciesNames), "species loaded."
			if not self.email:
				print "Please enter a valid email address to download sequence data from GenBank\n"
				self.email = raw_input("Email: ")
				Entrez.email = self.email
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
				print "\nDELETE MODE"
				print "\tSpID - irreversibly delete a species, e.g. '0'"
				print "\tgene - irreversibly delete an entire gene (brings up gene choice prompt)"
				print "\toutput - write out downloaded sequences in FASTA format"
				print "Other modes: 'reload', 'trim', 'replace', 'merge'. Hit enter to continue.\n"
			inputSeq = raw_input("DNA Editing (delete): ")
			if inputSeq:
				try:
					if int(inputSeq) in range(len(self.speciesNames)):
						del self.sequences[int(inputSeq)]
						del self.speciesNames[int(inputSeq)]
						if self.taxonomy:
							del self.taxonomy[int(inputSeq)]
						print "SeqID", inputSeq, "Successfully deleted"
						print "Re-calulating summary statistics..."
						self.dnaChecking()
						return 'delete', False
				except:
					pass
				if inputSeq == "gene":
					if len(self.genes) == 1:
						print "You only have one gene loaded!"
					else:
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
					return 'delete', False
				elif inputSeq == 'output':
					if self.sequences:
						for i,gene in enumerate(self.genes):
							currentGene = []
							if self.sequences[i]:
								for j,seq in enumerate(self.sequences[i]):
									if self.sequences[i][j]:
										currentGene.append(self.sequences[i][j])
							SeqIO.write(currentGene, self.stem+"_WIP_"+gene+".fasta", 'fasta')
					print "Raw sequences written out!"
					return 'delete', False
				elif inputSeq == "trim":
					return "trim", True
				elif inputSeq == "reload":
					return "reload", True
				elif inputSeq == "replace":
					return "replace", True
				elif inputSeq == "merge":
					return "merge", True
				else:
					print "Sorry,", inputSeq, "was not recognised. Please try again."
					return 'delete', False
			else:
				return "EXIT", True
		
		def reloadMode(firstTime=True):
			if not self.email:
				print "\nPlease enter a valid email address to download sequence data from Genbank\n"
				self.email = raw_input("Email: ")
				Entrez.email = self.email
			if firstTime:
				print "\nRELOAD MODE"
				print "\tSpID - reload all sequences for one species, e.g. '0'"
				print "\tSeqID GeneName - reload one gene for one species, e.g. '0rbcL'"
				print "\t>X / <X - reload all sequences longer/shorter than X, e.g. '>900' / '<900'"
				print "\tEVERYTHING - reload all sequences (search will be more thorough)"
				print "Other modes: 'delete', 'trim', 'replace', 'merge'. Hit enter to continue.\n"
			inputSeq = raw_input("DNA Editing (reload):")
			if inputSeq:
				try:
					if int(inputSeq) in range(len(self.speciesNames)):
						for i,each in enumerate(self.sequences[int(inputSeq)]):
							self.APICheck()
							self.sequences[int(inputSeq)][i] = sequenceDownload(self.speciesNames[int(inputSeq)], self.genes[i], thorough=True, retMax=self.genBankRetMax, seqChoice=self.replaceSeqChoice, targetLength=self.dnaCheck['quantileLengths'][i][2])
						print "Re-calulating summary statistics..."
						self.dnaChecking()
						return 'reload', False
				except:
					pass
				for i,gene in enumerate(self.genes):
					if gene in inputSeq:
						try:
							seqID = re.search("[0-9]*", inputSeq).group()
							seqID = int(seqID)
							if seqID and seqID < len(self.sequences):
								print "Reloading SeqID", seqID, "gene", gene
								self.APICheck()
								self.sequences[seqID][i] = sequenceDownload(self.speciesNames[seqID], self.genes[i], thorough=True, retMax=self.genBankRetMax, seqChoice=self.replaceSeqChoice, targetLength=self.dnaCheck['quantileLengths'][i][2])
								print "Re-calulating summary statistics..."
								self.dnaChecking()
								return 'reload', False
						except:
							pass
				if ">" in inputSeq:
					threshold = int(inputSeq.replace('>', ''))
					for i,sp in enumerate(self.sequences):
						for j,gene in enumerate(sp):
							if len(self.sequences[i][j]) > threshold:
								self.APICheck()
								self.sequences[i][j] = sequenceDownload(self.speciesNames[i], self.genes[j], thorough=True, retMax=self.genBankRetMax, seqChoice=self.replaceSeqChoice, targetLength=self.dnaCheck['quantileLengths'][j][2])
					print "Re-calulating summary statistics..."
					self.dnaChecking()
					return 'reload', False
				if "<" in inputSeq:
					threshold = int(inputSeq.replace('<', ''))
					for i,sp in enumerate(self.sequences):
						for j,gene in enumerate(sp):
							if len(self.sequences[i][j]) < threshold:
								self.APICheck()
								self.sequences[i][j] = sequenceDownload(self.speciesNames[i], self.genes[j], thorough=True, retMax=self.genBankRetMax, seqChoice=self.replaceSeqChoice, targetLength=self.dnaCheck['quantileLengths'][j][2])
					print "Re-calulating summary statistics..."
					self.dnaChecking()
					return 'reload', False
				if inputSeq == "EVERYTHING":
					geneOutput = findGenes(self.speciesNames, self.genes, seqChoice=self.replaceSeqChoice, verbose=True, thorough=True, download=True, retMax=self.genBankRetMax, targetNoGenes=self.nGenes, delay=self.delay, spacer=self.spacer)
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
				elif inputSeq == "merge":
					return "merge", True
				else:
					print "Sorry,", inputSeq, "was not recognised. Please try again."
					return 'reload', False
			else:
				return "EXIT", True
		
		def trimMode(firstTime=True):
			if firstTime:
				print "\nTRIM MODE"
				print"\tSpID - trim all that species' sequences, e.g., '0'"
				print "\tSpID GeneName - trim a single sequence for a single species, e.g., '0rbcL'"
				print "\t>X / <X - trims all sequences greater/lesser than X, e.g. '>1000' or '<1000'"
				print "\t'EVERYTHING' - trim all sequences"
				print "\t'type' - *IMPORTANT* select the type of gene (mitochondrial, nuclear, etc.) you've downloaded."
				print "First time you trim, I will use sequence annotations, the second time, I will search for ORFs. Check documentation for details "
				print "Other modes: 'delete', 'reload', 'replace', 'merge'. Hit enter to continue.\n"
			inputSeq = raw_input("DNA Editing (trim):")
			if inputSeq:
				try:
					if int(inputSeq) in range(len(self.speciesNames)):
						for i,each in enumerate(self.sequences[int(inputSeq)]):
							if self.sequences[int(inputSeq)]:
								self.sequences[int(inputSeq)][i] = cleanSequenceWrapper(self.sequences[int(inputSeq)][i], self.genes[i], DNAtype=self.codonModels[i], gapType='-')
						print "Re-calulating summary statistics..."
						self.dnaChecking()
						return 'trim', False
				except:
					pass
				for i,gene in enumerate(self.genes):
					if gene in inputSeq:
						try:
							seqID = re.search("[0-9]*", inputSeq).group()
							seqID = int(seqID)
							if seqID < len(self.sequences):
								print "Trimming SeqID", seqID, "gene", gene
								self.sequences[seqID][i] = cleanSequenceWrapper(self.sequences[int(inputSeq)][i], self.genes[i], DNAtype=self.codonModels[i], gapType='-')
								print "Re-calulating summary statistics..."
								self.dnaChecking()
								return 'trim', False
						except:
							pass
				if ">" in inputSeq:
					threshold = int(inputSeq.replace('>', ''))
					for i,sp in enumerate(self.sequences):
						for j,gene in enumerate(sp):
							if len(self.sequences[i][j]) > threshold:
								self.sequences[i][j] = cleanSequenceWrapper(self.sequences[i][j], self.genes[j], DNAtype=self.codonModels[j], gapType='-')
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
				elif inputSeq == "merge":
					return "merge", True
				elif inputSeq == "EVERYTHING":
					for i,sp in enumerate(self.sequences):
						for j,gene in enumerate(sp):
							if self.sequences[i][j]:
								self.sequences[i][j] = cleanSequenceWrapper(self.sequences[i][j], self.genes[j], DNAtype=self.codonModels[j], gapType='-', )
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
				print "\nPlease enter a valid email address to download sequence data from GenBank\n"
				self.email = raw_input("Email: ")
				Entrez.email = self.email
			if firstTime:
				print "\nREPLACE MODE"
				print "\tSpID - replace a species with a congener"
				print "\tEVERYTHING - replace all species without sequences with a congener"
				print "\tTHOROUGH - replace all species without sequences with a close-relative, where this doesn't conflict with a GenBank-derrived taxonomy of your species."
				print "Other modes: 'delete', 'reload', 'trim', 'merge'. Hit enter to continue.\n"
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
									self.speciesNames[i] = candidate[0] +" - " + self.speciesNames[i]
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
									replacements = cladeSpecies(lineage[0])
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
											self.speciesNames[i] = self.speciesNames[i]
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
									print "......alternative found:", candidate
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
								for candidate in self.uniqueTaxonomy[i][1:]:
									geneList = []
									locker = False
									for gene in self.genes:
										self.APICheck()
										temp = sequenceDownload(candidate, gene)
										geneList.append(temp)
										if temp:
											self.speciesNames[i] = self.speciesNames[i]
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
				elif inputSeq == "merge":
					return "merge", True
				else:
					print "Sorry,", inputSeq, "was not recognised. Please try again."
					return 'replace', False
			else:
				return "EXIT", True
		
		def mergeMode(firstTime=True):
			if firstTime:
				print "\nMERGE MODE"
				print "\tSpID,SpID - merge species into a group, e.g. 0,1,2"
				print "Groups must have only one species with sequence data. There is no 'undo' option, and you *must not* merge two merged groups."
				print "Other modes: 'delete', 'reload', 'trim', 'replace'. Hit enter to continue.\n"
			inputSeq = raw_input("DNA Editing (merge): ")
			if inputSeq:
				try:
					spp = [int(x) for x in inputSeq.split(',')]
					spp.sort()
					if len(spp) > 1:
						foundDNA = (False, -1)
						for i,sp in reversed(list(enumerate(spp))):
							if sp in range(len(self.sequences)):
								for j,gene in enumerate(self.sequences[sp]):
									if self.sequences[sp][j]:
										if foundDNA[0]:
											print "Sorry, you can't merge species if more than one has DNA. Please try again."
											return "merge", False
										else:
											foundDNA = (True, sp)
											break
							else:
								print "I couldn't find that SeqID."
								return 'merge', False
						if foundDNA[0]:
							#It's safe to do the deletions
							merged = [self.speciesNames[foundDNA[1]]]
							for i in reversed(spp):
								if i != foundDNA[1]:
									merged.append(self.speciesNames[i])
									del self.sequences[i]
									del self.speciesNames[i]
							self.mergedSpp.append(merged)
							print "Successfully merged species into", merged[0]
							print "Re-calulating summary statistics..."
							self.dnaChecking()			
							return "merge", False
						else:
							print "One of the species you're merging must have DNA data"
							return "merge", False
					else:
						print "You can't merge a single species! Maybe try 'delete' mode?"
						return "merge", False
				except:
					pass
				if inputSeq == "delete":
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
			elif mode == "merge":
				mode, firstTime = mergeMode(firstTime)
			else:
				raise RuntimeError("Unrecognised DNA Editing mode!")
			if mode:
				if mode == "EXIT":
					self.alignment = []
					locker = False
	
	def unmerge(self):
		def unmergeNewick(tree, sppList):
			Phylo.write(tree, 'mergedTempPhylogeny.tre', 'newick')
			with open('mergedTempPhylogeny.tre', 'r') as f:
				treeText = ''
				for each in f:
					treeText += each.strip()
			for i,spp in enumerate(sppList):
				sp = spp[0]
				length = float(re.search('(?<='+sp+'\:)[0-9\.]*', treeText).group())
				replacement = ':' + str(length/2) + ','
				insertion = '(' + replacement.join(spp)
				insertion += ':' + str(length/2) + '):' + str(length/2)
				treeText = re.sub('('+sp+'\:)[0-9\.]*', insertion, treeText)
			with open('mergedTempPhylogeny.tre', 'w') as f:
				f.write(treeText)
			mergedTree = Phylo.read('mergedTempPhylogeny.tre', 'newick')
			os.remove('mergedTempPhylogeny.tre')
			return mergedTree
		
		def unmergeBEAST(treeFile, speciesNames, sppList):
			with open(treeFile, 'r') as f:
				for each in f:
					if 'tree TREE1 = [&R]' in each:
						treeText = each.strip()
						treeText = re.sub('\[.*?\]','', treeText)
						treeText = treeText.replace('tree TREE1 = ', '')
						treeText = treeText.replace(' ', '')
			for i,sp in enumerate(speciesNames):
				treeText = re.sub('(?<=,|\()('+str(i+1)+')(?!([0-9]))', sp, treeText)
				with open(treeFile+"SMOOTH_TEMP", 'w') as f:
					f.write(treeText)
				tree = Phylo.read(treeFile+"SMOOTH_TEMP", 'newick')
				os.remove(treeFile+"SMOOTH_TEMP")
			return unmergeNewick(tree, sppList)
		
		
		if self.mergedSpp:
			if 'BEAST' in self.phylogenyMethods:
				self.phylogenyMerged = unmergeBEAST(self.phylogeny, self.speciesNames, self.mergedSpp)
				if self.smoothPhylogeny:
					self.smoothPhylogenyMerged = unmergeBEAST(self.smoothPhylogeny, self.mergedSpp)
				print "Warning: merged BEAST phylogenies do not have posteriors. For these, use the RAW output."
			else:
				self.phylogenyMerged = unmergeNewick(self.phylogeny, self.mergedSpp)
				if self.smoothPhylogeny:
					self.smoothPhylogenyMerged = unmergeNewick(self.smoothPhylogeny, self.mergedSpp)
			print "Smoothed tree successfully merged. Be aware that floating point precision issues may alter the tree's branchlengths."
		else:
			return False
	
	def alignmentEditing(self):
		alignmentDisplay(self.alignment, self.alignmentMethods, self.genes, checkAlignmentList(self.alignment, method='everything'))
		print "\n\t'output' - write out alignments. I recommend you look at your alignment before continuing"
		print "\t'DNA' - return to DNA editing stage"
		print "\t'align' - return to alignment stage, discarding current alignments."
		print "\t'trimal' - automatically trim your sequences using trimAl"
		print "\t'raxml=X' - run X RAxML runs for each alignment, and calculate the R-F distances between the trees and alignments"
		print "\t'metal' - calculate SSP distances between genes and alignments using metal"
		print "Hit enter to continue and choose one final alignment per gene\n"
		locker = True
		while locker:
			inputAlign = raw_input("Alignment Checking:")
			if inputAlign:
				if inputAlign == 'output':
					for i,gene in enumerate(self.genes):
						for j,method in enumerate(self.alignmentMethods):
							AlignIO.write(self.alignment[i][j], self.stem+"_"+gene+"_"+method+".fasta", 'fasta')
					print "...output written!"
				elif inputAlign == 'DNA':
					self.alignment = []
					self.alignmentMethod = False
					self.dnaChecking()
					self.dnaEditing()
					self.align()
					alignmentDisplay(self.alignment, self.alignmentMethods, self.genes, checkAlignmentList(self.alignment, method='everything'))
				elif inputAlign == 'align':
					self.alignmentMethod = False
					self.alignment = []
					self.align()
					alignmentDisplay(self.alignment, self.alignmentMethods, self.genes, checkAlignmentList(self.alignment, method='everything'))
				elif inputAlign == 'trimal':
					self.alignment = cleanAlignment(self.alignment, timeout=99999)
					alignmentDisplay(self.alignment, self.alignmentMethods, self.genes, checkAlignmentList(self.alignment, method='everything'))
				elif inputAlign == 'metal':
					if len(self.alignmentMethods) == 1 and len(self.genes) == 1:
						print "You've only used one alignment method, and one gene!"
					else:
						self.metal = metal(self.alignment)
						print "\nSSP distances between genes and alignments:"
						header = []
						for i,gene in enumerate(self.genes):
							if len(gene) > 4:
								gene = gene[0:4]
							else:
								gene = gene.ljust(4)
							for j,method in enumerate(self.alignmentMethods):
								header.append(gene + "_" + method[0:4])
						print " ".join(header)
						x = 1
						for i,gene in enumerate(self.metal):
							curRow = ["        "]
							x += 1
							for k,method in enumerate(gene):
								curRow.append(str(round(method, 4)).ljust(9))
							print " ".join(curRow)
				elif 'raxml' in inputAlign:
					if len(self.alignmentMethods) == 1 and len(self.genes) == 1:
						print "You've only used one alignment method, and one gene!"
					else:
						try:
							noTrees = int(inputAlign.replace('raxml=', ''))
							#Make starting trees and random seeds
							#startingTrees = []
							#for i in range(noTrees):
							#	startingTrees.append(RAxML(self.alignment[0][0], method='localVersion-startingOnly'))

							#Find trees from this search
							# - if you start with the same tree for each, it tends to be shit...
							trees = []
							for i,gene	in enumerate(self.alignment):
								for k,method in enumerate(gene):
									for tree in range(noTrees):
										trees.append(RAxML(method, 'localVersion'))

							#Calculate mean R-F distances
							#...between alignment methods
							Phylo.write(trees, 'alignCheckTemp.tre', 'newick')
							pipe = TerminationPipe('raxml -f r -z alignCheckTemp.tre -n alignCheckTemp -m GTRGAMMA', 999999)
							pipe.run()
							if not pipe.failure:
								self.alignRF = []
								with open('RAxML_RF-Distances.alignCheckTemp') as f:
									for line in f:
										temp = line.strip()
										self.alignRF.append(re.search("[0-9]{1}\.[0-9]+", temp).group())

								os.remove('alignCheckTemp.tre')
								os.remove('RAxML_RF-Distances.alignCheckTemp')
								os.remove('RAxML_info.alignCheckTemp')
								print "\nMean Robinson-Folds distances between genes and alignments:"
								header = []
								for i,gene in enumerate(self.genes):
									if len(gene) > 4:
										gene = gene[0:4]
									else:
										gene = gene.ljust(4)
									for j,method in enumerate(self.alignmentMethods):
										header.append(gene + "_" + method[0:4])
								header.insert(0, "         ")
								print " ".join(header)
								colToWrite = len(self.genes) * len(self.alignmentMethods) -1
								x = 0
								spacer = 1
								for row in range((len(self.genes) * len(self.alignmentMethods))-1):
									currRow = []
									for col in range(colToWrite):
										temp = []
										for rep in range(noTrees):
											temp.append(float(self.alignRF[x]))
											x += 1
										currRow.append(str(round(sum(temp)/len(temp), 3)).ljust(9))
									for i in range(spacer):
										currRow.insert(0, "         ")
									print header[row+1], " ".join(currRow)
									colToWrite -= 1
									spacer += 1
							else:
								print "!!!Something went wrong with the RAxML runs. Sorry!"
						except:
							print "Sorry, try again - something like 'raxml=5'."
				else:
					print "Sorry, I don't understand", inputAlign, "- please try again."
			else:
				locker = False
		if len(self.alignmentMethods) > 1:
			print "\nIn the prompt below is the name of a gene. Type the ID number of the alignment you'd like to use for that gene."
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
		else:
			for i,gene in enumerate(self.genes):
				self.alignment[i] = self.alignment[i][0]
	
	def align(self):
		methods = ['muscle', 'mafft', 'clustalo', 'prank', 'everything', 'quick']
		if not self.alignmentMethod:
			print "\nChoose one alignment method ('muscle', 'mafft', 'clustalo', 'prank'), or..."
			print "\t'everything' -  all four and compare their outputs"
			print "\t'quick' to do only the first three"
			print "Return will use MAFFT; prank is very slow!\n"
			locker = True
			while locker:
				alignInput = raw_input("DNA Alignment: ")
				if alignInput:
					if alignInput in methods:
						self.alignmentMethod = alignInput
						locker = False
					else:
						print "Sorry, I didn't recognise", alignInput, "- please try again."
				else:
					self.alignmentMethod = 'mafft'
					locker = False
		print "Starting alignment..."
		self.alignment = alignSequences(seqList=self.sequences, sppNames=self.speciesNames, method=self.alignmentMethod, tempStem='temp', timeout=99999999, nGenes=len(self.genes))
		print "\nAlignment complete!"
		if self.alignmentMethod == 'everything':
			self.alignmentMethods = ['muscle', 'mafft', 'clustalo', 'prank']
		elif self.alignmentMethod == 'quick':
			self.alignmentMethods = ['muscle', 'mafft', 'clustalo']
		else:
			self.alignmentMethods = [self.alignmentMethod]
	
	def phylogen(self, method="raxml-localVersion"):
		def raxmlSetup(options=False):
			def parseOptions(inputStr):
				methods = ''
				if 'raxml' == inputStr:
					methods += 'RAxML-defaults'
				if 'integratedBootstrap=' in inputStr:
					temp = inputStr.split('-')
					for each in temp:
						if 'integratedBootstrap' in each:
							methods += each + '-'
							break
				if 'restart=' in inputStr:
					temp = inputStr.split('-')
					for each in temp:
						if 'restart' in each:
							methods += each + '-'
							break
				if 'partitions' in inputStr:
					methods += 'partitions-'
				return methods
			
			if len(self.alignment) > 1:
				align, partitions = self.concatenateSequences()
			else:
				align, partitions = self.alignment, None
			if options:
				if not parseOptions(options):
					print "ERROR! Your RAxML options (", options, ") were not recognised. Please re-enter below."
					self.phylogenyMethods = ''
				else:
					self.phylogenyMethods = 'RAxML-' + parseOptions(options)
			
			if not self.phylogenyMethods:
				print "RAXML:"
				print "\t 'integratedBootstrap=X' - conduct X number of bootstraps and a thorough ML search in one run (!)"
				print "\t 'restart=X' - conduct X number of full ML searches (!)"
				print "\t 'partitions' - concatenate all genes into a single partition (not the default)"
				print "Specify multiple options with hyphens (e.g., 'restart=5-partitions'), but do not mix options marked with '(!)'"
				print "Hit enter to conduct one search\n"
				raxmlLock = True
				while raxmlLock:
					raxmlInput = raw_input("Phylogeny Building (RAxML): ")
					if raxmlInput:
						methods = parseOptions(raxmlInput)
						if not methods:
							print "Sorry, I don't understand", raxmlInput, "- please try again."
						raxmlLock = False
					else:
						self.phylogenyMethods = 'RAxML-'
						print "...running RAxML with default options:", self.phylogenyMethods
						raxmlLock = False
			
			self.phylogeny = RAxML(align, method=self.phylogenyMethods+'localVersion', constraint=self.constraint, timeout=999999, partitions=partitions)
		
		def beastSetup(options=False):
			def parseOptions(inputStr):
				methods = ''
				chainLength = 1000000
				logRate = 1000
				screenRate = 1000
				overwrite = True
				restart = 0
				burnin = 0.1
				if 'beast' == inputStr:
					methods += '-GTR-GAMMA'
				if 'GTR' in inputStr:
					methods +=	'-GTR'
				if 'HKY' in inputStr:
					methods += '-HKY'
				if 'GAMMA' in inputStr:
					methods +=	'-GAMMA'
				if 'overwriteBlock' in inputStr:
					overwrite = False
				#Oh my god make this simpler, it's absurd!!!
				if 'chainLength' in inputStr:
					temp = inputStr.split('-')
					for each in temp:
						if 'chainLength' in each:
							chainLength = int(each.replace('chainLength=', ''))
							methods += '-chainLength='+str(chainLength)
							break
				if 'logRate' in inputStr:
					temp = inputStr.split('-')
					for each in temp:
						if 'logRate' in each:
							logRate = int(each.replace('logRate=', ''))
							methods += '-logRate='+str(logRate)
							break
				if 'screenRate' in inputStr:
					temp = inputStr.split('-')
					for each in temp:
						if 'screenRate' in each:
							screenRate = int(each.replace('screenRate=', ''))
							methods += '-screenRate='+str(screenRate)
							break
				if 'restart' in inputStr:
					temp = inputStr.split('-')
					for each in temp:
						if 'restart' in each:
							restart = int(each.replace('restart=', ''))
							methods += '-restart='+str(restart)
							break
				if 'burnin' in inputStr:
					temp = inputStr.split('-')
					for each in temp:
						if 'burnin' in each:
							burnin = int(each.replace('burnin=', ''))
							methods += '-burnin='+str(burnin)
							break
				return methods, logRate, screenRate, chainLength, overwrite, burnin
			
			if options:
				if not parseOptions(options)[0]:
					print "ERROR! Your BEAST options (", options, ") were not recognised. Please re-enter below."
				else:
					self.phylogenyMethods, logRate, screenRate, chainLength, overwrite, burnin = parseOptions(options)
					self.phylogenyMethods = 'BEAST' + self.phylogenyMethods
			
			if not self.phylogenyMethods:
				print "BEAST:"
				print "\t 'GTR' - use a GTR model (default) (!)"
				print "\t 'HKY' - use an HKY model (!)"
				print "\t 'GAMMA' - use four gamma rate categories"
				print "\t 'chainLength=X' - use chain length of 'X' (default 1 000 000)"
				print "\t 'logRate=X' - log output every 'X' steps (default 1000)"
				#print "\t 'screenRate=X' - report ouptut to the screen every 'X' steps"
				print "\t 'burnin=X' - discard 'X' initial fraction of chain as a 'burnin' (default 0.1 = 10%)"
				print "\t 'overwriteBlock' - causes BEAST to halt if there are any files from previous attempts in your working directory"
				#print "\t 'restart=X' - conduct X independent searches"
				print "You must specify a DNA model (GTR or HKY) if defining your own parameters"
				print "Specify multiple options with hyphens (e.g., 'restart=5-partitions'), but do not mix options marked with '(!)'"
				print "Hit enter to use the defaults\n"
				beastLock = True
				while beastLock:
					beastInput = raw_input("Phylogeny Building (BEAST): ")
					methods = ''
					if beastInput:
						self.phylogenyMethods, logRate, screenRate, chainLength, overwrite, burnin = parseOptions(beastInput)
						self.phylogenyMethods = 'BEAST' + self.phylogenyMethods
						if self.phylogenyMethods:
							beastLock = False
						else:
							print "Sorry, I don't understand", beastInput, "- please try again."
					else:
						print "...running BEAST with default options"
						self.phylogenyMethods, logRate, screenRate, chainLength, overwrite, burnin = parseOptions('')
						self.phylogenyMethods = 'GTR-GAMMA'
						beastLock = False
			
			if not 'GTR' in self.phylogenyMethods and not 'HKY' in self.phylogenyMethods:
				self.phylogenyMethods += '-GTR-GAMMA'
			print "...running BEAST with options ", self.phylogenyMethods
			if len(self.alignment) > 1:
				self.phylogeny, self.beastXML, self.beastTrees, self.beastLogs = BEAST(self.alignment, method=self.phylogenyMethods, constraint=self.constraint, timeout=999999, chainLength=chainLength, logRate=logRate, screenRate=screenRate, overwrite=overwrite, burnin=burnin)
			else:
				self.phylogeny, self.beastXML, self.beastTrees, self.beastLogs  = BEAST(self.alignment[0], method=self.phylogenyMethods, constraint=self.constraint, timeout=999999, chainLength=chainLength, logRate=logRate, screenRate=screenRate, overwrite=overwrite, burnin=burnin)
		
		if self.phylogenyMethods:
			if 'beast' in self.phylogenyMethods:
				beastSetup(self.phylogenyMethods)
			elif 'raxml' in self.phylogenyMethods:
				raxmlSetup(self.phylogenyMethods)
			else:
				print "I don't understand your phylogeny construction method", self.phylogenyMethods, ". Please enter one now."
				self.phylogenyMethods = ''
		
		if not self.phylogenyMethods:
			print "You can either build a maximum likelihood tree ('raxml') or a Bayesian tree ('beast'). If unsure, hit enter to use RAxML\n"
			locker = True
			while locker:
				phyloInput = raw_input("Phylogeny Building: ")
				if phyloInput:
					if phyloInput == 'raxml':
						raxmlSetup()
						locker = False
					elif phyloInput == 'beast':
						beastSetup()
						locker = False
					else:
						print "Sorry, I don't understand", phyloInput, "- please try again."
				else:
					print "...using RAxML..."
					raxmlSetup('')
					locker= False
	
	def rateSmooth(self):
		def PATHd8():
			if sys.platform == 'win32':
				try:
					subprocess.Popen('/requires/PATHd8')
				except:
					print "\n*ERROR*\n"
					print "PATHd8 requires 'cygwin' to be installed, and you don't have it."
					print "'website' - open cygwin website and install it"
					print "...or... anything else to return to the rate-smoothing prompt"
					cygwin = raw_input('Rate-smoothing (Cygwin): ')
					if cygwin and cygwin == 'website':
						browser.open('http://www.cygwin.com/')
					else:
						return False
			print "\nPlease enter an outgroup for your phylogeny, 'species' to see the tips in your phylogeny, or just hit enter to cancel and continue."
			pathd8Locker = True
			spNames = [x.name for x in self.phylogeny.get_terminals()]
			while pathd8Locker:
				inputPathd8 = raw_input("Rate-Smoothing (PATHd8): ")
				if inputPathd8:
					if inputPathd8 == 'species':
						print "\nSpecies names:"
						for each in spNames:
							print each
					elif inputPathd8 in spNames:
						self.phylogeny.root_with_outgroup(inputPathd8)
						length = 0
						if type(self.alignment) is list:
							for i,align in enumerate(self.alignment):
								length += align.get_alignment_length()
						else:
							length = self.alignment.get_alignment_length()
						self.smoothPhylogeny = rateSmooth(self.phylogeny, sequenceLength=length)
						pathd8Locker = False
						print "...phylogeny rate-smoothed! Continuing..."
					else:
						print "Sorry, I couldn't find", inputSmooth, " in your phylogeny - try again"
				else:
					print "...cancelling rate-smoothing and continuing..."
					pathd8Locker = False
			return True
		
		def smoothBEAST():
			print "You're about to perform a BEAST search with the topology constrained to that of your best phylogen(y/ies)."
			print "This means you have a choice of BEAST options:"
			print "\t 'GTR' - conduct search with a GTR model (!)"
			print "\t 'HKY' - conduct search with an HKY model (!)"
			print "\t 'GAMMA' - conduct search with four gamma rate categories"
			print "\t 'chainLength=X' - conduct search with chain length 'X'"
			print "\t 'logRate=X' - log output every 'X' steps"
			#print "\t 'screenRate=X' - report ouptut to the screen every 'X' steps"
			print "\t 'burnin=X' - Discard 'X' initial fraction of chain as a 'burnin', e.g. 0.1 = 10%"
			print "\t 'overwriteBlock' - causes BEAST to halt if there are any files from previous attempts in your working directory"
			print "...you *must* specify a DNA model (at least GTR or HKY) if defining your own parameters"
			print "...to specify multiple options, type them all separated by hyphens (e.g. 'GTR-GAMMA')"
			print "...the options with '(!)' after them cannot be used in conjunction with each other"
			print "...or... just hit enter to use the defaults!"
			beastLock = True
			while beastLock:
				beastInput = raw_input("Rate-Smoothing (BEAST): ")
				methods = ''
				chainLength = 1000000
				logRate = 1000
				screenRate = 1000
				overwrite = True
				burnin = 0.1
				if beastInput:
					if 'GTR' in beastInput:
						methods +=	'-GTR'
					if 'HKY' in beastInput:
						methods += '-HKY'
					if 'GAMMA' in beastInput:
						methods +=	'-GAMMA'
					if 'overwriteBlock' in beastInput:
						overwrite = False
					if 'chainLength' in beastInput:
						temp = beastInput.split('-')
						for each in temp:
							if 'chainLength' in each:
								chainLength = int(each.replace('chainLength=', ''))
								methods += '-chainLength='+str(chainLength)
								break
					if 'logRate' in beastInput:
						temp = beastInput.split('-')
						for each in temp:
							if 'logRate' in each:
								logRate = int(each.replace('logRate=', ''))
								methods += '-logRate='+str(logRate)
								break
					if 'screenRate' in beastInput:
						temp = beastInput.split('-')
						for each in temp:
							if 'screenRate' in each:
								screenRate = int(each.replace('screenRate=', ''))
								methods += '-screenRate='+str(screenRate)
								break
					if 'burnin' in beastInput:
						temp = beastInput.split('-')
						for each in temp:
							if 'burnin' in each:
								burnin = int(each.replace('burnin=', ''))
								methods += '-burnin='+str(burnin)
								break
					if methods:
						self.rateSmoothMethods = 'BEAST' + methods
						if len(self.alignment) > 1:
							self.smoothPhylogeny, self.smoothBeastXML, self.smoothBeastTrees, self.smoothBeastLogs  = BEAST(self.alignment, method=self.rateSmoothMethods, constraint=self.phylogeny, logRate=logRate, screenRate=screenRate, chainLength=chainLength, overwrite=overwrite, burnin=burnin, timeout=999999, tempStem='beast_smooth')
						else:
							self.smoothPhylogeny, self.smoothBeastXML, self.smoothBeastTrees, self.smoothBeastLogs = BEAST(self.alignment[0], method=self.rateSmoothMethods, constraint=self.phylogeny, logRate=logRate, screenRate=screenRate, chainLength=chainLength, overwrite=overwrite, timeout=999999, burnin=burnin, tempStem='beast_smooth')
						beastLock = False
					else:
						print "Sorry, I don't understand", beastInput, "- please try again."
				else:
					print "...running BEAST with default options"
					self.rateSmoothMethods = 'BEAST-GTR-GAMMA'
					if len(self.alignment) > 1:
						self.smoothPhylogeny, self.smoothBeastXML, self.smoothBeastTrees, self.smoothBeastLogs = BEAST(self.alignment, method=self.rateSmoothMethods, constraint=self.phylogeny, logRate=logRate, screenRate=screenRate, chainLength=chainLength, overwrite=overwrite, timeout=999999, burnin=burnin, tempStem='beast_smooth')
					else:
						self.smoothPhylogeny, self.smoothBeastXML, self.smoothBeastTrees, self.smoothBeastLogs = BEAST(self.alignment[0], method=self.rateSmoothMethods, constraint=self.phylogeny, logRate=logRate, screenRate=screenRate, chainLength=chainLength, overwrite=overwrite, timeout=999999, burnin=burnin, tempStem='beast_smooth')
					beastLock = False
					
		
		print "\nYou can now rate-smooth your phylogeny (i.e., make its branch lengths proportional to evolutionary time)."
		print "\t'pathd8' - rate-smooth using PATHd8 (fastest!)"
		print "\t'beast' - rate-smooth using BEAST (doesn't require you to specify an outgroup)"
		print "...or hit enter to continue without rate smoothing."
		locker = True
		while locker:
			inputSmooth = raw_input("Rate-Smoothing: ")
			if inputSmooth:
				if inputSmooth == 'pathd8':
					success = PATHd8()
					if success:
						locker = False
				elif inputSmooth == 'beast':
					smoothBEAST()
					locker = False
				else:
					print "Sorry, I can't ", inputSmooth, " - try again."
			else:
				print "...continuing without rate-smoothing"
				locker = False
	
	def cleanUpSequences(self):
		cleaned = []
		for i,sp in reversed(list(enumerate(self.sequences))):
			foundSequence = False
			for j,seq in reversed(list(enumerate(sp))):
				if seq:
					foundSequence = True
				else:
					self.sequences[i][j] = SeqRecord(Seq(""))
					self.sequences[i][j].id = self.speciesNames[i] + "Empty" + str(i)
					self.sequences[i][j].name = self.speciesNames[i] + "Empty" + str(i)
					self.sequences[i][j].description = 'Empty sequence made up by phyloGenerator'
			
			if not foundSequence:
				cleaned.append(self.speciesNames[i])
				del self.sequences[i]
				del self.speciesNames[i]
				if self.taxonomy:
					del self.taxonomy[i]
		
		if cleaned:
			print "\nThe following species did not have any DNA associated with them, and so have been excluded:"
			for each in cleaned:
				print each
	
	def renameSequences(self):
		for i,name in enumerate(self.speciesNames):
			self.speciesNames[i] = name.replace(" ", "_")
		
		if self.mergedSpp:
			for i,sppList in enumerate(self.mergedSpp):
				for j,sp in enumerate(sppList):
					self.mergedSpp[i][j] = sp.replace(" ", "_")
		
		for i in range(len(self.sequences)):
			tGenBankIDs = []
			for k in range(len(self.genes)):
				if self.sequences[i][k]:
					tGenBankIDs.append(self.sequences[i][k].id)
					self.sequences[i][k].id = self.speciesNames[i]
				else:
					tGenBankIDs.append("NO_SEQUENCE")
			self.genBankIDs.append(tGenBankIDs)
	
	def writeOutput(self):
		#Change working directory
		os.chdir(self.workingDirectory)
		#Log - TO-DO
		#Sequences
		if self.sequences:
			for i,gene in enumerate(self.genes):
				currentGene = []
				for seq in self.sequences:
					currentGene.append(seq[i])
				SeqIO.write(currentGene, self.stem+"_"+gene+".fasta", 'fasta')
		
		#Alignment
		if type(self.alignment) is list:
			for i,align in enumerate(self.alignment):
				AlignIO.write(align, self.stem+"_"+self.genes[i]+"_alignment.fasta", 'fasta')
		else:
			AlignIO.write(self.alignment, self.stem+"_alignment.fasta", 'fasta')
		
		#Sequence info
		if self.genBankIDs:
			for i,gene in enumerate(self.genes):
				with open(self.stem+"_"+self.genes[i]+"_sequence_info.txt", 'w') as f:
					f.write("Species Name, Sequence ID\n")
					for j,name in enumerate(self.speciesNames):
						f.write(name + "_" + self.genBankIDs[j][i] + "\n")
		
		#Phylogeny
		if self.phylogenyMerged:
			if 'BEAST' in self.phylogenyMethods:
				os.rename(self.oldDirectory + '/' + self.phylogeny, self.stem+"_"+self.genes[i]+"_RAW_phylogeny.nex")
				os.rename(self.oldDirectory + '/' + self.beastXML, self.stem+"_"+self.genes[i]+"_RAW_BEAST.xml")
				os.rename(self.oldDirectory + '/' + self.beastTrees, self.stem+"_"+self.genes[i]+"_RAW_BEAST.trees")
				os.rename(self.oldDirectory + '/' + self.beastLogs, self.stem+"_"+self.genes[i]+"_RAW_BEAST.log")
				Phylo.write(self.phylogenyMerged, self.stem+"_MERGED_phylogeny.tre", 'newick')
			else:
				Phylo.write(self.phylogenyMerged, self.stem+"_MERGED_phylogeny.tre", 'newick')
				Phylo.write(self.phylogeny, self.stem+"_RAW_phylogeny.tre", 'newick')
			
		else:
			if self.phylogeny:
				if 'BEAST' in self.phylogenyMethods:
					os.rename(self.oldDirectory + '/' + self.phylogeny, self.stem+"_"+self.genes[i]+"_phylogeny.nex")
					os.rename(self.oldDirectory + '/' + self.beastXML, self.stem+"_"+self.genes[i]+"_RAW_BEAST.xml")
					os.rename(self.oldDirectory + '/' + self.beastTrees, self.stem+"_"+self.genes[i]+"_RAW_BEAST.trees")
					os.rename(self.oldDirectory + '/' + self.beastLogs, self.stem+"_"+self.genes[i]+"_RAW_BEAST.log")
				else:
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
						f.write(",".join(each) + "\n")
		
		#Smoothed phylogeny
		if self.smoothPhylogeny:
			if self.smoothPhylogenyMerged:
				if 'BEAST' in self.rateSmoothMethods:
					os.rename(self.oldDirectory + '/' + self.smoothPhylogeny, self.stem+"_"+self.genes[i]+"_RAW_smoothed_phylogeny.nex")
					os.rename(self.oldDirectory + '/' + self.smoothBeastXML, self.stem+"_"+self.genes[i]+"_RAW_smoothed_BEAST.xml")
					os.rename(self.oldDirectory + '/' + self.smoothBeastTrees, self.stem+"_"+self.genes[i]+"_RAW_smoothed_BEAST.trees")
					os.rename(self.oldDirectory + '/' + self.smoothBeastLogs, self.stem+"_"+self.genes[i]+"_RAW_smoothed_BEAST.log")
					Phylo.write(self.smoothPhylogenyMerged, self.stem+"_"+self.genes[i]+"_MERGED_smoothed_phylogeny.tre", 'newick')
				else:
					Phylo.write(self.smoothPhylogeny, self.stem+"_RAW_smoothed_phylogeny.tre", 'newick')
					Phylo.write(self.smoothPhylogenyMerged, self.stem+"_MERGED_smoothed_phylogeny.tre", 'newick')
			else:
				if 'BEAST' in self.rateSmoothMethods:
					os.rename(self.oldDirectory + '/' + self.smoothPhylogeny, self.stem+"_"+self.genes[i]+"_RAW_smoothed_phylogeny.nex")
					os.rename(self.oldDirectory + '/' + self.smoothBeastXML, self.stem+"_"+self.genes[i]+"_RAW_smoothed_BEAST.xml")
					os.rename(self.oldDirectory + '/' + self.smoothBeastTrees, self.stem+"_"+self.genes[i]+"_RAW_smoothed_BEAST.trees")
					os.rename(self.oldDirectory + '/' + self.smoothBeastLogs, self.stem+"_"+self.genes[i]+"_RAW_smoothed_BEAST.log")
				else:
					Phylo.write(self.smoothPhylogeny, self.stem+"_smoothed_phylogeny.tre", 'newick')
	
	def getConstraint(self, fileName=""):
		def newickConstraint():
			if self.constraintFile:
				try:
					self.constraint = Phylo.read(self.constraintFile, 'newick')
				except IOError:
					print "\nConstraint tree not found. Error!"
					sys.exit()
			else:
				print "\nEnter the location of your *unrooted* constraint tree (in newick format)"
				print "Press enter to cancel\n"
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
							print "Constraint tree does *not* match the species names you've inputed - this is case sensitive. Please load another file."
					else:
						print "...No constraint tree loaded"
						locker = False
						return False
		
		def phylomatic():
			if self.phylomaticTaxonomy and self.phylomaticPhylogeny:
				try:
					self.constraint = unrootPhylomatic(createConstraintTree(self.phylomaticTaxonomy, method="phylomaticTaxonomy", fileName=self.phylomaticPhylogeny, tempStem='temp'))
					return False
				except:
					print "Phylomatic didn't work - check your input files"
					print "Returning you to constraint method choice prompt...\n"
					self.phylomaticTaxonomy, self.phylomaticPhylogeny = False, False
			
			print "\nCreating a constraint tree using Phylomatic."
			if not self.phylomaticPhylogeny:
				print "Enter the filename of the reference phylogeny (in newick format) below. Hit enter to abort"
				locker = True
				while locker:
					inputPhylogeny = raw_input("Phylomatic (reference tree): ")
					if inputPhylogeny:
						try:
							open(inputPhylogeny)
							print ""
							self.phylomaticPhylogeny = inputPhylogeny
							locker = False
						except IOError:
							print "\nFile not found. Please try again!"
					else:
						print "...No constraint tree loaded. Continuing."
						return False
			if not self.phylomaticTaxonomy:
				print "Enter your species' taxonomy (in Phylomatic's format)\n"
				locker = True
				while locker:
					inputTaxonomy = raw_input("Phylomatic (taxonomy): ")
					if inputTaxonomy:
						try:
							open(inputTaxonomy)
							self.phylomaticTaxonomy = inputTaxonomy
							locker = False
						except IOError:
							print "\nFile not found. Please try again!"
					else:
						print "...No constraint tree loaded. Continuing."
						locker = False
			try:
				self.constraint = unrootPhylomatic(createConstraintTree(self.phylomaticTaxonomy, method="phylomaticTaxonomy", fileName=self.phylomaticPhylogeny, tempStem='temp'))
				return False
			except:
				print "Returning you to constraint method choice prompt...\n"
				return True
		
		def taxonomy():
			if not self.email:
				print "\nPlease enter a valid email address to download sequence data from GenBank\n"
				self.email = raw_input("Email: ")
				Entrez.email = self.email
			print "\nCreating a 'taxonomy' for your species from GenBank"
			print "\tNote: this will be saved as a file, but will not generate a constraint tree (yet...)"
			for sp in self.speciesNames:
				self.APICheck()
				self.taxonomy.append(findLineage(sp))
			print "...lineages found!"
		
		stopper = True
		if not self.constraintMethod:
			print "I recommend you use a constraint tree with this program"
			print "\t'newick' - supply your own constraint tree"
			print "\t'phylomatic' - use Phylomatic to generate a tree"
			print "\t'taxonomy' - download the NCBI taxonomy for your species (does not generate a constraint tree)"
			print "Warning: Phylomatic can trim the end off species names, causing conflicts with phyloGenerator that are hard to detect. Rooted phylogenies are *not* valid constraint."
			print "Otherwise, press enter to continue without a constraint tree.\n"
			stopper = True
			while stopper:
				constraintInput = raw_input("Constraint Method: ")
				if constraintInput:
					if constraintInput == 'newick':
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
					return False
		else:
			if self.constraintMethod == 'newick':
				newickConstraint()
			elif self.constraintMethod == 'phylomatic':
				phylomatic()
			elif self.constraintMethod == 'taxonomy':
				taxonomy()
			else:
				print "Error in constraint tree method name."
				#do something better here please!
				return False
		if self.constraint:
			try:
				if len(self.genes) > 1:
					align, partitions = self.concatenateSequences()
				else:
					align = self.alignment
					partitions = self.partitions
				RAxML(align, constraint=self.constraint)
			except:
				os.remove('RAxML_info.tempOut')
				print "\nYour constraint phylogeny is incompatible with RAxML, which probably means you used a rooted phylogeny with Phylomatic."
				print "Your constraint has been written out as 'temp_contraint.tre'. You can edit it if you wish."
				return True
			print "To conduct RAxML searches with and without your constraint tree, and calculate the Robinson-Foulds distances between them, enter the number of times you would like to do a tree search below. Otherwise, simply press enter."
			checkLocker = True
			while checkLocker:
				checkerInput = raw_input("Constraint Method (check): ")
				if checkerInput:
					try:
						nSearch = int(checkerInput)
						trees = []
						for i in range(nSearch):
							trees.append(RAxML(align, partitions=partitions))
						for i in range(nSearch):
							trees.append(RAxML(align, partitions=partitions, constraint=self.constraint))
						Phylo.write(trees, 'constraintCheckTemp.tre', 'newick')
						pipe = TerminationPipe('raxml -f r -z constraintCheckTemp.tre -n constraintCheckTemp -m GTRGAMMA', 999999)
						pipe.run()
						if not pipe.failure:
							with open('RAxML_RF-Distances.constraintCheckTemp') as f:
								for line in f:
									temp = line.strip()
									self.constraintRFs.append(float(re.search("[0-9]{1}\.[0-9]+", temp).group()))
							os.remove('constraintCheckTemp.tre')
							os.remove('RAxML_RF-Distances.constraintCheckTemp')
							os.remove('RAxML_info.constraintCheckTemp')
							aFree = False
							freeFree = []
							freeConstrained = []
							constrainedConstrained = []
							x = 0
							for i in range((nSearch*2) -1):
								aFree = not aFree
								bFree = not aFree
								for j in range(i, (nSearch*2) -1):
									if aFree:
										if bFree:
											freeFree.append(self.constraintRFs[x])
										else:
											freeConstrained.append(self.constraintRFs[x])
									else:
										if bFree:
											freeConstrained.append(self.constraintRFs[x])
										else:
											constrainedConstrained.append(self.constraintRFs[x])
									x += 1
									bFree = not bFree
							print "\tConstrained mean distance:   ", str(round(mean(constrainedConstrained),2)), " (SD: ", str(round(std(constrainedConstrained), 4)), ")"
							print "\tUnconstrained mean distance: ", str(round(mean(freeFree),2)), " (SD: ", str(round(std(freeFree), 4)), ")"
							print "\tMean distance between them:  ", str(round(mean(freeConstrained),2)), " (SD: ", str(round(std(freeConstrained), 4)), ")"
						checkLocker = False
					except:
						print "Sorry, I didn't get that. Please try again."
				else:
					checkLocker = False
	
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
			partitions = [0, self.alignment[0].get_alignment_length()]
			tempAlignment = self.alignment[0]
			for i in range(1, (len(self.genes)-1)):
				tempAlignment += self.alignment[i]
				partitions.append(lengths[-1] + self.alignment[i])
			tempAlignment
			for i,x in enumerate(self.speciesNames):
				tempAlignment[i].id = self.speciesNames[i].replace(' ', '_')
		return tempAlignment, partitions
	
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
		print "\nDNA CHECKING"
		currentState.dnaChecking()
		print "You may now edit the sequences you are using. Deleting species may change species' IDs"
		print "Huge variation in lengths of sequences (e.g., thousands of base pairs) crashes many alignment programs"
		print "All species without sequence data will be ignored when continuing to the next step."
		currentState.dnaEditing()
		
		#DNA Cleanup and renaming
		currentState.cleanUpSequences()
		currentState.renameSequences()
		
		#Alignment
		print "\nDNA ALIGNMENT"
		currentState.align()
		print "\nALIGNMENT CHECKING"
		currentState.alignmentEditing()
		
		#Constraint tree
		print "\nCONSTRAINT TREE"
		currentState.getConstraint()
		
		#Phylogeny building
		print "\nPHYLOGENY BUILDING"
		currentState.phylogen()
		
		#Rate smoothing
		if 'BEAST' in currentState.phylogenyMethods:
			print "\nSKIPPING RATE SMOOTHING STEP"
			print "\t(unecessary with BEAST phylogeny)"
		else:
			print "\nRATE SMOOTHING"
			currentState.rateSmooth()
		
		#Handle merged species
		currentState.unmerge()
		
		#Output
		currentState.writeOutput()
		print "\nCongratulations! Exiting phyloGenerator."

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="phyloGenerator - phylogeny generation for ecologists.", epilog="Help at http://willpearse.github.com/phyloGenerator - written by Will Pearse")
	parser.add_argument("--version", action="store_true", help="Display version information.")
	parser.add_argument("--manual", action="store_true", help="(Attempt to) open browser and show help")
	parser.add_argument("-name", "-n", help="'Stem' name for all output files.")
	parser.add_argument("-wd", help="Working directory for all output files")
	parser.add_argument("-dna", "-d", help="Unaligned DNA (in FASTA format).")
	parser.add_argument("-alignment", "-a", help="Alignment method")
	parser.add_argument("-phylogen", "-p", help="Phylogeny construction method and options")
	parser.add_argument("-gene", "-g", help="The genes to search for (multiple genes are comma-separated)")
	parser.add_argument("-nGenes", "-ng", help="The number of genes to search for (if fewer than suggested in 'genes' are required)")
	parser.add_argument("-species", "-s", help="Binomial names of species, each on a new line")
	parser.add_argument("-consTree", "-cT", help="Constraint tree (in newick format).")
	parser.add_argument("-consPhylomat", "-cP", help="Phylomatic's phylogeny (comma) and taxonomy")
	parser.add_argument("-email", "-e", help="Email address for GenBank searches.")
	parser.add_argument("-options", "-o", help="Options file giving detailed instructions to phyloGen.")
	parser.add_argument("-delay", help="Delay (seconds) when pausing between any internet API calls.")
	main()
