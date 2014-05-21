#!/usr/bin/env python                                                                  

'''                                                                                     
djc 4/22/14                                                                         
'''

import sys
import csv

def convert(predsFileName, outFileName):
    	"""Convert from my pred file format to gff"""

    	with open(predsFileName) as reads, open(outFileName, 'w') as outFile:
        	predsReader = csv.reader(reads, delimiter='\t')
		writer = csv.writer(outFile, delimiter='\t')

		chrome = ""
		lastEnd = 0
		first = True

		for row in predsReader:	
			if(first):
				first = False
				continue

			if(row[0] != chrome):
				lastEnd = 0
				chrome = row[0]

			out = ("CHROMOSOME_"+row[0], "Custom", "Prediction", lastEnd, row[1], ".", ".", ".", "Prediction="+row[3]+";Number of sites="+row[4]+";N2/Genotyped="+row[5]+";MY14/Genotyped="+row[6]+";Heterozygous/Genotyped="+row[7]+";Unknown/All="+row[8])

			writer.writerow(out)

			lastEnd = row[2]
		

if __name__ == '__main__':
	convert(*sys.argv[1:])
	
