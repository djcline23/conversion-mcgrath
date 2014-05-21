#!/usr/bin/env python                                                                  

'''                                                                                     
djc 3/31/14                                                                         
'''

import sys
import csv
import scipy.stats
from collections import defaultdict

OLD_PRINT = False

class Region:
	prevPred = "-"
	nextPred = "-"
	
	def __init__(self, strain, chrome, start, earlyEnd, lateEnd, pred, sites):
		self.strain = strain
		self.chrome = chrome
		self.start = start
		self.earlyEnd = earlyEnd
		self.lateEnd = lateEnd
		self.pred = pred
		
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
		if OLD_PRINT:
			return '\t'.join((self.chrome, str(self.earlyEnd), str(self.lateEnd), self.pred, str(self.numSites), str(self.percentN), str(self.percentM), str(self.percentH), str(self.percentU)))
		else:
			return '\t'.join((self.strain, self.chrome, "{:,}".format(self.start), "{:,}".format(self.earlyEnd), "{:,}".format(self.size), "{:,}".format(self.numSites), "{:.2%}".format(self.percentN), "{:.2%}".format(self.percentM), "{:.2%}".format(self.percentH), "{:.2%}".format(self.percentI), "{:.2%}".format(self.percentU), self.prevPred, self.pred, self.nextPred, self.link))

def cutoffs(binSize):
	"""Figure out cutoffs for the given bin size.
	Returns a tuple of cutoffs. Less than the first number
	means the bin is most likely M, greater than the second
	number means the bin is most likely N, in between is 
	heterozygous"""
	
	binSize = int(binSize)	
	#Binomial for number of Ns in bin
	mProb = scipy.stats.binom.pmf(range(binSize+1),binSize,0.01)
	hProb = scipy.stats.binom.pmf(range(binSize+1),binSize,0.5)
	nProb = scipy.stats.binom.pmf(range(binSize+1),binSize,0.99)
	
	mThresh = -1
	nThresh = binSize+1

	for i in range(binSize+1):
		if(hProb[i] > mProb[i]):
			mThresh = i
			break

	for i in range(binSize, 0, -1):
		if(hProb[i] > nProb[i]):
			nThresh = i
			break;

	return float(mThresh), float(nThresh)

def predict(inFileName, outFileName, strain, binSize):
	"""Walks through a file of site calls and predicts
	Whether the region is M, N, or H based on a sliding
	window of the given binSize"""

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
		if OLD_PRINT:
			predFile.write('\t'.join(("Chromosome", "Breakpoint Lower", "Breakpoint Upper", "Prediction", "Number of SNPs", "N2/Genotyped", "MY14/Genotyped", "Heterozygous/Genotyped", "Unknown/All")) + '\n')

		else:
			predFile.write(','.join(("Strain", "Chromosome", "Start", "End", "Size (kbp)", "Number of SNPs", "N2/Genotyped", "MY14/Genotyped", "Heterozygous/Genotyped", "Incongruities/Genotyped", "Unknown/All", "Previous Prediction", "Prediction", "Next Prediction", "Link")) + '\n')

		#Gather up thresholds for the binSize and binSize/2 because sometimes
		#we break the window in half for more positional knowledge
		binSize = int(binSize)
		fullBinThresh = cutoffs(binSize)
		halfBinThresh = cutoffs(binSize / 2)
		thresholds = (fullBinThresh[0], fullBinThresh[1], halfBinThresh[0], halfBinThresh[1])

		#Make predictions for each chromesome
		for chrome in sorted(chromeMap.keys()):
			predictChrome(strain, chrome, chromeMap[chrome], binSize, thresholds, predFile)

def predictChrome(strain, chrome, sites, binSize, thresholds, predFile):
	"""Walks through the sites and makes predictions"""

	currPred = None
	currStart = 0
	currEnd = binSize

	#Keep track of where the last cut point was to make summary stats when
	#we close the current prediction
	lastCut = 0
	
	#Link previous and next predictions
	prevRegion = None
	
	while currEnd < len(sites):
		window = sites[currStart:currEnd]

		newPred, halfCall = call(window, thresholds, currPred)

		if currPred and currPred != newPred:
			#Then we have had a transition

			if(halfCall):
				#Only look in first half of window for cutoff
				window = window[:(binSize/2)]

			#Find best cutoff point
			cut = findBestCut(window, currPred, newPred)

			currStart = currStart + cut

			region = Region(strain, chrome, getStart(lastCut, sites), sites[currStart - 1][1], sites[currStart][1], currPred, sites[lastCut:currStart - 1])
			processNewRegion(predFile, prevRegion, region)

			lastCut = currStart
			prevRegion = region

		currPred = newPred

		#Not considering overlapping windows, not sure if that's the right thing to do or not
		#sometimes we only call half of the window length, in which case only
		#advance half a window size
		if(halfCall):
			currEnd = currEnd - binSize / 2

		currStart = findNextDiff(sites, currEnd+1, currPred)
		currEnd = currStart + binSize

	#Now we are at the end, print the last prediction
	region = Region(strain, chrome, getStart(lastCut, sites), sites[-1][1], sites[-1][1], currPred, sites[lastCut:currStart - 1])
	processNewRegion(predFile, prevRegion, region)
	predFile.write(region.toString() + '\n')

def processNewRegion(outfile, prevRegion, newRegion):
	"""Sets previous/next prediction for the two regions and prints the previous region"""
	if(prevRegion):
		prevRegion.nextPred = newRegion.pred
		newRegion.prevPred = prevRegion.pred
		
		outfile.write(prevRegion.toString() + '\n')

def getStart(lastCut, sites):
	"""Gets the start of a region based on the last cutting point"""
	if(lastCut == 0):
		return 0
	else:
		return sites[lastCut][1]

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
	that the first half was N, and the seocnd half was M (or vice versa).
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
	#cutoffs(*sys.argv[1:])
	predict(*sys.argv[1:])
	
