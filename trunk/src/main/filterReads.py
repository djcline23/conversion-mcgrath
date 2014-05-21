#!/usr/bin/env python                                                                  

'''                                                                                     
djc 4/7/14                                                                         
'''

import sys
import csv

def filterReads(readsFileName, filterFileName, outFileName):
    	"""Only keeps rows from the given file of reads which have a position
	found in filter"""

    	with open(readsFileName) as reads, open(filterFileName) as filterFile, open(outFileName, 'w') as outFile:
        	readsReader = csv.reader(reads, delimiter='\t')
		filterReader = csv.reader(filterFile, delimiter='\t')
		writer = csv.writer(outFile, delimiter='\t')

		filterLine = next(filterReader, None)

		for row in readsReader:	
			#If we hit the end of the filter file we are done
			if not filterLine:
				break;

			if( (row[0] == filterLine[0]) and (row[1] == filterLine[1]) ):
				writer.writerow(row)
				filterLine = next(filterReader, None)
		

if __name__ == '__main__':
	filterReads(*sys.argv[1:])
	
