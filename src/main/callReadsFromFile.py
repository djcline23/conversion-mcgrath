#!/usr/bin/env python                                                                  

'''                                                                                     
djc 3/27/14                                                                         
'''

import sys
import csv
import re
import util

#Call het if at least one of each if True. Call het based on binomials if False
SIMPLE_HET_CALL = True

def callReads(cutoffs, mutationsFileName, pileupFileName, outFileName):
  """For each position in the mutations file, run mpileup and                         
    parse to make a call for that position across all reads. A read                     
    can be called as N2, MY14, or pool quality. For now poor quality                    
    just means a base that does not match N2 or MY14.                                   
                                                                                        
    After calling all the reads, we make a call for the position as                     
    N2 if it only has N2 reads, MY14 if it only has MY14 reads,                         
    Heterozygous if it has at least one read of each type, or                           
    Unknown if there are no reads or only poor quality reads."""
    
  with open(mutationsFileName) as mutations, open(pileupFileName) as pileup, open(outFileName, 'w') as outFile:
    mutReader = csv.reader(mutations, delimiter='\t')
    pileReader = csv.reader(pileup, delimiter='\t')

    currPileLine = next(pileReader, None)

    for row in mutReader:    
      # parse and call each site
      site = 'U'
      n = 0
      m = 0
      p = 0
  
      siteChrom  = row[3]
      sitePos    = row[4]
  
      # If there is no pileup line, no need for checking
      if currPileLine:
        siteN      = row[5]
        siteM      = row[6]
        pileChrom  = currPileLine[0].split("_")[1]
        pilePos    = currPileLine[1]
        pileStr    = currPileLine[4]
  
        # Some sites might not appear in the pileup file.
        # Because both files are ordered, if the current line from the pileup
        # is not for the same location as the current site, we can just call it
        # U with no reads and go to the next site
        # print siteChrom
        # print pileChrom
        # print sitePos
        # print pilePos
  
        if(siteChrom == pileChrom and sitePos == pilePos):    
          reads = clean(pileStr.upper())
          n = reads.count(siteN)
          m = reads.count(siteM)
          p = len(reads) - n - m
          currPileLine = next(pileReader, None)
      
        site = call(n, m, cutoffs)

    outFile.write('\t'.join((siteChrom, sitePos, site, str(n), str(m), str(p))) + '\n')

def call(n, m, cutoffs):
  call = 'U'
  
  if SIMPLE_HET_CALL:
    if(n > 0 and m > 0):
      call = 'H'
    elif(n > 0):
      call = 'N'
    elif(m > 0):
      call = 'M'
  else:  
    total = n + m
  
    if(total > 0):
      # The total number of reads determine how many instances we need to make a call
      thresholds = cutoffs[total]
    
      if(n > thresholds[1]):
        call = 'N'
      elif(n < thresholds[0]):
        call = 'M'
      else:
        call = 'H'
      
  return call

def clean(reads):
  """Clean up mpileup calls so that we are left with just the bases,
  not the extra stuff it adds"""

  # Easy stuff first
  reads = re.sub(r"[$*<>]|\^.", "", reads)

  # Now do the trickier inserts/deletes
  # First find the regex
  while(1):
    match = re.search(r"[+-](\d+)", reads)
    if not match:
      break

    num = int(match.group(1))
    # now find the index and keep the string before and after the insert/delete
    startIndex = reads.find(match.group(0))
    endIndex = startIndex + 1 + len(match.group(1)) + num
    reads = reads[0:startIndex] + reads[endIndex:]    

  return reads

if __name__ == '__main__':
  # print clean('*A$C>G<T+23ACACACACACACACACACAGGGGT^JC-2TTC-11AAAAAAAAAAAC+4TTTTC')
  hetCutoffs = {}
  
  for i in range(1, 500):
    hetCutoffs[i] = util.cutoffs(i)
    
  callReads(hetCutoffs, *sys.argv[1:])
  
