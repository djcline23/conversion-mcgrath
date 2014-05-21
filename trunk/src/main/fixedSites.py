#!/usr/bin/env python                                                                  

'''                                                                                     
djc 4/7/14                                                                         
'''

import os
import sys
import csv
from collections import OrderedDict

def findFixed(sitesFileName, readsDir):
	"""For each position in the sites file, check whether it is
	found to have been called as both MY14 and N2 within some
	set of reads in the reads directory. If so, save it to notFixed.txt,
	else save to fixed.txt"""

	#For each site, keep track of whether we have seen a read for both MY14 and N2
	candSites = OrderedDict()
	notFixed = []

    	with open(sitesFileName) as sites:
        	siteReader = csv.reader(sites, delimiter='\t')

		for line in siteReader:
			#tuple of chrome, pos as key for dict, then whether we have seen M
			#and N respectively
			#Use ints for position so sorting works properly for notFixed
			candSites[(line[0].split("_")[1], int(line[1]))] = (False, False)

		#With that dictionary, now check the sites in each file
		for fileName in os.listdir(readsDir):
			print fileName
			with open(os.path.join(readsDir, fileName)) as currFile:
				reader = csv.reader(currFile, delimiter='\t')

				for line in reader:
					key = (line[0], int(line[1]))

					#If we have removed this site from consideration, 
					#keep going
					if not key in candSites:
						continue

					call = line[2]					
					check = True

					if(call == 'M'):
						candSites[key] = (True, candSites[key][1])
					elif(call == 'N'):
						candSites[key] = (candSites[key][0], True)
					elif(call == 'H'):
						candSites[key] = (True, True)
					else:
						check = False

					#If that site has seen both now, move it to notFixed
					if check and candSites[key][0] and candSites[key][1]:
						del candSites[key]
						notFixed.append(key)

						#If we found both for every site, quit looking
						if(len(candSites) == 0):
							break

		print "Fixed: " + str(len(candSites))
		print "Not Fixed: " + str(len(notFixed))

	with open('fixed.txt', 'w') as fixedFile, open('notFixed.txt', 'w') as notFixedFile:
		fixedWriter = csv.writer(fixedFile, delimiter='\t')
		for site in candSites.keys():
			fixedWriter.writerow(site)

		#Sort by position first
		notFixed.sort()
		notFixedWriter = csv.writer(notFixedFile, delimiter='\t')	
		for site in notFixed:
			notFixedWriter.writerow(site)

if __name__ == '__main__':
	#print clean('*A$C>G<T+23ACACACACACACACACACAGGGGT^JC-2TTC-11AAAAAAAAAAAC+4TTTTC')
	findFixed(*sys.argv[1:])
	
