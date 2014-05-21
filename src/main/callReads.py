#!/usr/bin/env python                                                                  

'''                                                                                     
djc 3/27/14                                                                         
'''

import sys
import csv
import subprocess
import re

def callReads(mutationsFileName, bamFileName, outFileName):
    """For each position in the mutations file, run mpileup and                         
    parse to make a call for that position across all reads. A read                     
    can be called as N2, MY14, or pool quality. For now poor quality                    
    just means a base that does not match N2 or MY14.                                   
                                                                                        
    After calling all the reads, we make a call for the position as                     
    N2 if it only has N2 reads, MY14 if it only has MY14 reads,                         
    Heterozygous if it has at least one read of each type, or                           
    Unknown if there are no reads or only poor quality reads."""

    with open(mutationsFileName) as mutations, open(outFileName, 'w') as outFile:
        mutReader = csv.reader(mutations, delimiter='\t')
        
	for row in mutReader:
		#do pileup
		pileup = subprocess.check_output(['samtools', 'mpileup', '-r', 'CHROMOSOME_%s:%s-%s' % (row[3], row[4], row[4]), bamFileName], stderr=subprocess.PIPE) 
		
		#parse and call
		site = 'U'
		n = 0
		m = 0
		p = 0
		if(len(pileup) is not 0):		
			reads = clean(pileup.split('\t')[4].upper())
			n = reads.count(row[5])
			m = reads.count(row[6])
			p = len(reads) - n - m

		if(n > 0 and m > 0):
			site = 'H'
		elif(n > 0):
			site = 'N'
		elif(m > 0):
			site = 'M'

		outFile.write('\t'.join((row[3], row[4], site, str(n), str(m), str(p))) + '\n')

def clean(reads):
	"""Clean up mpileup calls so that we are left with just the bases,
	not the extra stuff it adds"""

	#Easy stuff first
	reads = re.sub(r"[$*<>]|\^.", "", reads)

	#Now do the trickier inserts/deletes
	#First find the regex
	while(1):
		match = re.search(r"[+-](\d+)", reads)
		if not match:
			break

		num = int(match.group(1))
		#now find the index and keep the string before and after the insert/delete
		startIndex = reads.find(match.group(0))
		endIndex = startIndex + 1 + len(match.group(1)) + num
		reads = reads[0:startIndex] + reads[endIndex:]		

	return reads

if __name__ == '__main__':
	#print clean('*A$C>G<T+23ACACACACACACACACACAGGGGT^JC-2TTC-11AAAAAAAAAAAC+4TTTTC')
	callReads(*sys.argv[1:])
	
