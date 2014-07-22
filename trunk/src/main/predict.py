#!/usr/bin/env python                                                                  

'''                                                                                     
djc 3/31/14                                                                         
'''

import math
import argparse
import csv
import util
from collections import defaultdict

#Table format (true) or gff format (false)
TABLE_PRINT = False

#Whether to delete chromosomes which contain a het region
DELETE_HET = False	

#Whether to extends region end points to the midpoint between genotyped sites or not
EXTEND_ENDS = False

DEFAULT_BINSIZE = 100

chromeLengths = {"I" : 15072000, "II" : 15279000, "III" : 13784000, "IV" : 17494000, "V" : 20920000, "X" : 17719000}

class Region:
	prevPred = "-"
	nextPred = "-"
	
	def __init__(self, strain, chrome, prevRegion, start, earlyEnd, lateEnd, pred, sites):
		self.strain = strain
		self.chrome = chrome
		self.start = start
		self.earlyEnd = earlyEnd
		self.lateEnd = lateEnd
		self.pred = pred
		
		if(prevRegion):
			prevRegion.nextPred = pred
			self.prevPred = prevRegion.pred
		
		#Calculate summary stats
		self.size = (earlyEnd - start) / 1000
		
		self.numSites = 0
		n = 0
		m = 0
		h = 0
		u = 0
	
		for site in sites:
			#Number of Us, plus this site
			self.numSites = self.numSites + site[2] + 1
	
			u = u + site[2]
	
			if(site[0] == 'N'):
				n = n + 1
			elif(site[0] == 'M'):
				m = m+1
			else:
				h = h + 1
	
		#The U counts are for sites before the given n/m/h site
		#So don't count for the first site unless it is the start
		#of the chromosome
		if start != 0:
			self.numSites = self.numSites - sites[0][2]
			u = u - sites[0][2]
	
		#For N,M,H do percentage of genotyped
		self.percentN = n / float(self.numSites - u)
		self.percentM = m / float(self.numSites - u)
		self.percentH = h / float(self.numSites - u)
		self.percentU = u / float(self.numSites)
		
		#Incongruities
		self.percentI = 0.0
		
		if(pred == 'M'):
			self.percentI = self.percentN + self.percentH
		elif(pred == 'N'):
			self.percentI = self.percentM + self.percentH
			
		#Show some context with the link
		#Less than 0 is bad, but greater than chrome size is okay
		#size is in kilobases, so convert to bases = *1000, then take a tenth = *100
		#at least 200 bp
		padding = max(self.size * 100, 1000)
		displayStart = max(self.start - padding, 0)
		displayEnd = self.earlyEnd + padding
			
		self.link = "<a target=\"_blank\" href=\"http://gbrowse2014.biology.gatech.edu/jbrowse.html?data=conversion&loc=CHROMOSOME_{0}%3A{1}..{2}&tracks=DNA%2C{3}_calledReads.H%2C{3}_calledReads.M%2C{3}_calledReads.N%2C{3}_calledReads.U%2C{3}_pred_100\">View</a>".format(self.chrome, str(displayStart), str(displayEnd), self.strain)
			
	def toString(self):
		if TABLE_PRINT:
			return '\t'.join((self.chrome, str(self.earlyEnd), str(self.lateEnd), self.pred, str(self.numSites), str(self.percentN), str(self.percentM), str(self.percentH), str(self.percentU)))
		else:
			return '\t'.join((self.strain, self.chrome, "{:,}".format(self.start), "{:,}".format(self.earlyEnd), "{:,}".format(self.size), "{:,}".format(self.numSites), "{:.2%}".format(self.percentN), "{:.2%}".format(self.percentM), "{:.2%}".format(self.percentH), "{:.2%}".format(self.percentI), "{:.2%}".format(self.percentU), self.prevPred, self.pred, self.nextPred, self.link))

def predict(inFileName, outFileName, strain, binSize):
	"""Walks through a file of site calls and predicts
	Whether the region is M, N, or H based on a sliding
	window of the given binSize.
	Returns a tuple with the number of chromosomes and regions predicted"""

	with open(inFileName) as sitesFile, open(outFileName, 'w') as predFile:
		#Load whole file into memory since it is only a couple of MBs
		reader = csv.reader(sitesFile, delimiter='\t')
		#dictionary of lists
		chromeMap = defaultdict(list)

		#We really only want to work with sites that have been called as N,M,H
		#so just count up the number of U sites between the more meaningful
		#sites for summary stats later
		uCount = 0
		for row in reader:
			if(row[2] == 'U'):
				uCount = uCount + 1
			else:
				chromeMap[row[0]].append((row[2], int(row[1]), uCount)) 
				uCount = 0

		#Write header
		if TABLE_PRINT:
			predFile.write(','.join(("Strain", "Chromosome", "Start", "End", "Size (kbp)", "Number of SNPs", "N2/Genotyped", "MY14/Genotyped", "Heterozygous/Genotyped", "Incongruities/Genotyped", "Unknown/All", "Previous Prediction", "Prediction", "Next Prediction", "Link")) + '\n')
		else:
			predFile.write('\t'.join(("Chromosome", "Breakpoint Lower", "Breakpoint Upper", "Prediction", "Number of SNPs", "N2/Genotyped", "MY14/Genotyped", "Heterozygous/Genotyped", "Unknown/All")) + '\n')
			
		#Gather up thresholds for the binSize and binSize/2 because sometimes
		#we break the window in half for more positional knowledge
		fullBinThresh = util.cutoffs(binSize)
		halfBinThresh = util.cutoffs(binSize / 2)
		thresholds = (fullBinThresh[0], fullBinThresh[1], halfBinThresh[0], halfBinThresh[1])

		#The total number of chromosomes and regions predicted
		chromeCount = 0
		regCount = 0

		#Make predictions for each chromosome
		for chrome in sorted(chromeMap.keys()):
			chromes, regs = predictChrome(strain, chrome, chromeMap[chrome], binSize, thresholds, predFile)
			chromeCount += chromes
			regCount += regs
			
		return chromeCount, regCount

def predictChrome(strain, chrome, sites, binSize, thresholds, predFile):
	"""Walks through the sites and makes predictions, printing them out to predFile.
	Also returns a tuple with:
	 0 or 1 in the first position for whether this chromosome was included in our predictions
	 the number of regions predicted on this chromosome"""

	currPred = None
	currStart = 0
	currEnd = binSize
	regions = []
	
	#If we find a het region, toss the whole chromosome
	foundHet = False

	#Keep track of where the last cut point was to make summary stats when
	#we close the current prediction
	lastCut = 0
	
	#Link previous and next predictions
	prevRegion = None
	
	while currEnd < len(sites):
		window = sites[currStart:currEnd]

		newPred, halfCall = call(window, thresholds, currPred)
		
		if(DELETE_HET and newPred == 'H'):
			foundHet = True
			break

		if currPred and currPred != newPred:
			#Then we have had a transition

			if(halfCall):
				#Only look in first half of window for cutoff
				window = window[:(binSize/2)]

			#Find best cutoff point
			cut = findBestCut(window, currPred, newPred)

			currStart = currStart + cut

			currRegion = Region(strain, chrome, prevRegion, getStart(lastCut, sites), getEarlyEnd(currStart, sites), getLateEnd(currStart, sites), currPred, sites[lastCut:currStart - 1])
			regions.append(currRegion)

			lastCut = currStart
			prevRegion = currRegion

		currPred = newPred

		#Advance the minimum of spaces that are guaranteed to be the same as the call - thresholds[3]
		nextSearchStart = max(currStart + int(math.ceil(thresholds[3])), lastCut)

		currStart = findNextDiff(sites, nextSearchStart, currPred)
		currEnd = currStart + binSize

	if(DELETE_HET and foundHet):
		print "Heterozygous region on {} chromosome {}, ignoring entire chromosome.".format(strain, chrome)
		return 0, 0
	else:
		#Now we are at the end, add the last prediction
		regions.append(Region(strain, chrome, prevRegion, getStart(lastCut, sites), chromeLengths[chrome], chromeLengths[chrome], currPred, sites[lastCut:currStart - 1]))
		
		for region in regions:
			predFile.write(region.toString() + '\n')
			
		return 1, len(regions)
	
	
def midpoint(position, sites):
	return (sites[position - 1][1] + sites[position][1]) / 2

def getStart(start, sites):
	"""Gets the start of a region based on the last cutting point"""
	if(start == 0):
		return 0
	else:
		if EXTEND_ENDS:
			return midpoint(start, sites)
		else:
			return sites[start][1]
		
def getEarlyEnd(end, sites):
	if EXTEND_ENDS:
		return midpoint(end, sites)
	else:
		return sites[end - 1][1]
	
def getLateEnd(end, sites):
	if EXTEND_ENDS:
		return midpoint(end, sites)
	else:
		return sites[end][1]

def findNextDiff(sites, index, curr):
	"""Returns the index of the next site that is different from the given
	one, or an index out of range if no different site is found"""
	while index < len(sites):
		if sites[index][0] != curr:
			break
		else:
			index = index + 1

	return index

def findBestCut(sequence, oldPred, newPred):
	"""Finds the cut point to transition from oldPred to newPred"""
	bestCut = 0
	bestWrong = len(sequence) + 1

	for i in range(len(sequence)):
		wrong = callWrong(sequence[:i], oldPred) + callWrong(sequence[i:], newPred)

		if wrong < bestWrong:
			bestCut = i
			bestWrong = wrong
			
	return bestCut

def callWrong(sequence, pred):
	"""Returns the number of sites that disagree with the given prediction"""

	wrong = 0

	#No wrong calls for hetero
	if(pred != 'H'):
		for site in sequence:
			if(site[0] != pred):
				wrong = wrong + 1

	return wrong


def call(window, thresholds, currPred):
	ret = callWindow(window, thresholds[0], thresholds[1])

	if(ret == 'H'):
		#Look more closely into het calls
		return callHetero(window, thresholds, currPred)
	else:
		return ret, False
		

def callWindow(window, mThresh, nThresh):
	nVotes = 0.0

	for site in window:
		if site[0] == 'N':
			nVotes = nVotes + 1
		elif site[0] == 'H':
			nVotes = nVotes + 0.5

	if nVotes < mThresh:
		return 'M'
	elif nVotes > nThresh: 
		return 'N'
	else:
		return 'H'

def callHetero(window, thresholds, currPred):
	"""We were finding that het calls were too liberal in the case
	that the first half was N, and the second half was M (or vice versa).
	So now divide a het window in half and proceed based on the content of
	each half."""

	half = len(window) / 2;
	first = window[:half]

	firstCall = callWindow(first, thresholds[2], thresholds[3])

	if(firstCall == 'H'):
		#What we call depends on the second window
		second = window[half:]
		secondCall = callWindow(second, thresholds[2], thresholds[3])

		if(secondCall == 'H'):
			#Okay we believe you, the whole window is het
			return 'H', False
		elif(secondCall == currPred):
			#This looks like XHX, so call the first part as truly het
			return 'H', True
		else:
			#This looks like XHY, so assume this is really just a transition
			#from X to Y. Call the whole thing as Y and let other code
			#determine the transition point
			return secondCall, False
	else:
		#Then take the call, but only for this first part
		return firstCall, True

if __name__ == '__main__':	
	parser = argparse.ArgumentParser()
	parser.add_argument("strain", help="the strain found in callsFile")
	parser.add_argument("callsFile", help="the file with site calls")
	parser.add_argument("outputFile", help="the file to save predictions to")
	parser.add_argument("-d",  action="store_true", help="delete heterozygous chromosomes")
	parser.add_argument("-e",  action="store_true", help="extend regions to touch each other")
	parser.add_argument("-t",  action="store_true", help="print in the format needed for making the table, as opposed to the gffs")
	parser.add_argument("-b",  "--binSize", type=int, help="bin size for predictions")
	args = parser.parse_args()
	
	DELETE_HET = args.d
	EXTEND_ENDS = args.e
	TABLE_PRINT = args.t
	
	binSize = DEFAULT_BINSIZE
	if args.binSize:
		binSize = args.binSize
	
	predict(args.callsFile, args.outputFile, args.strain, binSize)
	
