#!/usr/bin/env python                                                                  

'''                                                                                     
djc 7/15/14                                                                         
'''

import predict

import os
import sys

if __name__ == '__main__':
  inFolder, outFolder, binSize = sys.argv[1:]
  
  chromeCount = 0
  regCount = 0
  
  for readFile in os.listdir(inFolder):
    if readFile.endswith(".txt"):
      base = readFile.replace("_calledReads.filtered.txt", "")
      outFile = os.path.join(outFolder, base + ".pred." + binSize + ".txt")
      readFile = os.path.join(inFolder, readFile)
      #print outFile
      
      chromes, regs = predict.predict(readFile, outFile, base, binSize)
      chromeCount += chromes
      regCount += regs
      
  print "Chromosomes: " + str(chromeCount)
  print "Regions: " + str(regCount)
