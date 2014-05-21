#!/usr/bin/env python                                                                  

'''                                                                                     
djc 4/7/14                                                                         
'''

import sys
import csv

def filterReads(pileFileName, filterFileName, outFileName):
    	"""Only keeps rows from the given pileup file which have a position
	found in filter"""

    	with open(pileFileName) as piles, open(filterFileName) as filterFile, open(outFileName, 'w') as outFile:
        	pileReader = csv.reader(piles, delimiter='\t')
		filterReader = csv.reader(filterFile, delimiter='\t')
		writer = csv.writer(outFile, delimiter='\t')

		good = {(loc[0], loc[1]) for loc in filterReader}

		for row in pileReader:	
			chrom = row[0].split("_")[1]

			if( (chrom, row[1]) in good ):
				row[0] = chrom
				writer.writerow(row)
		

if __name__ == '__main__':
	filterReads(*sys.argv[1:])
	
