#!/usr/bin/env python                                                                  

'''                                                                                     
djc 3/31/14                                                                         
'''

import sys
import csv
import scipy.stats
from collections import defaultdict

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

	print mThresh
	print nThresh

	return float(mThresh), float(nThresh)

def predict(inFileName, outFileName, binSize):
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
				chromeMap[row[0]].append((row[2], row[1], uCount)) 
				uCount = 0

		#Write header
		predFile.write('\t'.join(("Chromosome", "Breakpoint Lower", "Breakpoint Upper", "Prediction", "Number of SNPs", "N2/Genotyped", "MY14/Genotyped", "Heterozygous/Genotyped", "Unknown/All")) + '\n')

		#Make predictions for each chromesome
		for chrome in sorted(chromeMap.keys()):
			predictChrome(chrome, chromeMap[chrome], int(binSize), predFile)

def predictChrome(chrome, sites, binSize, predFile):
	"""Walks through the sites and makes predictions"""

	#Thresholds for predicting based on frequency of calls in the window
	mThresh, nThresh = cutoffs(binSize)

	currPred = None
	currStart = 0
	currEnd = binSize

	#Keep track of where the last cut point was to make summary stats when
	#we close the current prediction
	lastCut = 0
	
	while currEnd < len(sites):
		window = sites[currStart:currEnd]

		newPred = call(window, mThresh, nThresh)

		if currPred and currPred != newPred:
			#Then we have had a transition

			#Find best cutoff point
			cut = findBestCut(window, currPred, newPred)

			currStart = currStart + cut


			predFile.write('\t'.join((chrome, sites[currStart - 1][1], sites[currStart][1], currPred)) + '\t')
			predFile.write('\t'.join(summaryStats(sites[lastCut:currStart - 1], lastCut == 0)) + '\n')
			lastCut = currStart

		currPred = newPred

		#Not considering overlapping windows, not sure if that's the right thing to do or not
		currStart = findNextDiff(sites, currEnd+1, currPred)
		currEnd = currStart + binSize

	#Now we are at the end, print the last prediction
	predFile.write('\t'.join((chrome, sites[-1][1], sites[-1][1], currPred)) + '\t')
	predFile.write('\t'.join(summaryStats(sites[lastCut:currStart - 1], lastCut == 0)) + '\n')

def summaryStats(sites, isChromeStart):
	"""Returns a tuple with the number of sites and percentage of N,M,H,U"""
	total = 0
	n = 0
	m = 0
	h = 0
	u = 0

	for site in sites:
		#Number of Us, plus this site
		total = total + site[2] + 1

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
	if not isChromeStart:
		total = total - sites[0][2]
		u = u - sites[0][2]

	#For N,M,H do percentage of genotyped
	pn = n / float(total - u)
	pm = m / float(total - u)
	ph = h / float(total - u)
	pu = u / float(total)

	return (str(total), str(pn), str(pm), str(ph), str(pu))

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


def call(window, mThresh, nThresh):
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

if __name__ == '__main__':
	#cutoffs(*sys.argv[1:])
	predict(*sys.argv[1:])
	
