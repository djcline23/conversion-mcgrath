#!/usr/bin/env python                                                                  

'''                                                                                     
djc 4/24/14                                                                         
'''

import math
import os
import sys
import csv
from collections import defaultdict

#Share these arrays to store our stats
outliers = []
expSeqErrors = []
expBarCodeErrors = []

def gather(bcRateFileName, outFilename):
	"""Gathers the data needed for my stats project!
	# of outliers
	expected # of sequencing errors
	expected # of barcode errors"""

	with open(bcRateFileName) as bcFile, open(outFilename, 'w') as outFile:

		#we needs the pileup data, the calls, and the predictions
		callDir = "calls"
		pileDir = "piles"
		predDir = "preds"

		callFilenames = iter(os.listdir(callDir))
		pileFilenames = iter(os.listdir(pileDir))
		predFilenames = iter(os.listdir(predDir))

		#The entries in bcFile and all the call, pileup, pred files should be
		#synchronized

		maxCalls = 0.0
		#So each step here represents one sample
		for row in bcFile:
			currCall = os.path.join(callDir, callFilenames.next())
			currPile = os.path.join(pileDir, pileFilenames.next())
			currPred = os.path.join(predDir, predFilenames.next())

			numCalls = handleSample(currCall, currPile, currPred, float(row))

			if(numCalls > maxCalls):
				maxCalls = numCalls

		#print maxCalls
		#Now we have all stats, write em out!
		outFile.write('list(n={},\noutliers=c({}),\nexpSeqErrors=c({}),\nexpBarcodeErrors=c({}))'.format(len(outliers), ", ".join(outliers), ", ".join(expSeqErrors), ", ".join(expBarCodeErrors)))

def handleSample(callFilename, pileFilename, predFilename, bcRate):
	"""Gathers the stats for the given sample"""

	print callFilename

	with open(callFilename) as callFile, open(pileFilename) as pileFile, open(predFilename) as predFile:
		#load each file into memory and sort by chromosome
        	callChromes = defaultdict(list)
		pileChromes = defaultdict(list)
		predChromes = defaultdict(list)

		callReader = csv.reader(callFile, delimiter='\t')
		pileReader = csv.reader(pileFile, delimiter='\t')
		predReader = csv.reader(predFile, delimiter='\t')

		#Only care about the call itself
		for row in callReader:
			callChromes[row[0]].append((int(row[1]), row[2])) 

		#Only care about quality scores
		for row in pileReader:
			pileChromes[row[0]].append((int(row[1]), row[5]))

		#Need first 4 rows
		#Skip first row
		first = True
		for row in predReader:
			if first:
				first = False
				continue

			predChromes[row[0]].append((int(row[1]), int(row[2]), row[3]))

		#Now handle each chromosome
		maxCalls = 0.0
		for chrome in predChromes.keys():
			numCalls = handleChromosome(callChromes[chrome], pileChromes[chrome], predChromes[chrome], bcRate)
			if(numCalls > maxCalls):
				maxCalls = numCalls

		return maxCalls

def handleChromosome(calls, piles, preds, bcRate):
	"""Walks through prediction file and gathers stats for each region"""

	#We will need to walk through calls and pileups just once
	callPos = 0
	callLength = len(calls)
	pileLength = len(piles)
	pilePos = 0	

	currCall = calls[0]
	currPile = piles[0]

	#The chance of not having a barcode error
	chanceNoBc = 1.0 - bcRate

	regionStart = 0

	for region in preds:
		pred = region[2]
		if(pred == "H"): #don't care about het, so keep going
			continue;

		outlierCount = 0.0
		expBarcode   = 0.0
		expSeqCount  = 0.0

		regionEnd   = region[0]

		#print "\n"
		#print start
		#print end

		#Catch the iters up to the current region
		while(currCall[0] < regionStart):
			callPos = callPos + 1

			if(callPos < callLength):
				currCall = calls[callPos]
			else:
				break

		while(currPile[0] < regionStart):
			pilePos = pilePos + 1

			if(pilePos < pileLength):
				currPile = piles[pilePos]
			else:
				break

		#Now go through each and gather the stats
		numCalls = 0.0
		while(currCall[0] <= regionEnd and numCalls < 1000.0):
			#print currCall[0]
			if(currCall[1] != "U"):
				numCalls = numCalls + 1.0

			if(currCall[1] != "U" and currCall[1] != pred):
				outlierCount = outlierCount + 1
			
			callPos = callPos + 1

			if(callPos < callLength):				
				currCall = calls[callPos]
			else:
				break

		#print numCalls

		sites = 0

		while(currPile[0] <= regionEnd and sites < 1000):
			sites = sites + 1
			#chance that the base was sequenced incorrectly, divide by 3 that it was opposite base
			qualScores = currPile[1]
			expSeqCount = expSeqCount + (chanceWrong(qualScores) / 3.0)

			#The chance of having a barcode error is 1 minus the chance that all reads did not have a barcode error
			#Divide by 2 because it might be a barcode error but still have the same base
			chanceBc = (1.0 - pow(chanceNoBc, len(qualScores))) / 2.0
			#print chanceNoBc
			#print len(qualScores)
			#print chanceBc
			expBarcode = expBarcode + chanceBc

			pilePos = pilePos + 1

			if(pilePos < pileLength):
				currPile = piles[pilePos]
			else:
				break

		#Set start for next runthrough
		regionStart = region[1]

		#print numCalls

		#Only taking regions that are at least 1000 data points
		if(outlierCount < 100.0 and numCalls > 999 and sites > 999):
			outliers.append(str(outlierCount).upper())
			expSeqErrors.append(str(expSeqCount).upper())
			expBarCodeErrors.append(str(expBarcode).upper())

		return numCalls

def chanceWrong(qualStr):
	"""Returns the chance that this base was sequence incorrectly based on phred scores of reads"""
	
	chanceRight = 1.0
	
	for char in qualStr:
		score = ord(char) - 33
		chanceRight = chanceRight * (1.0 - math.pow(10.0, score / -10.0))

	return 1.0 - chanceRight

if __name__ == '__main__':
	gather(*sys.argv[1:])
	
