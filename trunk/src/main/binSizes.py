#!/usr/bin/env python                                                                  

'''                                                                                     
djc 7/11/14                                                                         
'''

import csv

def binRegions(filename):
  bins = [0 for i in range(500)]
  
  with open(filename) as regionsFile:
    reader = csv.reader(regionsFile, delimiter='\t')
    
    for row in reader:
      binRegions = int(row[4].replace(",", "")) / 50
      bins[binRegions] = bins[binRegions] + 1
      
    #Bin 419 corresponds to longest chromosome
    for i in range(419):
      print bins[i]

if __name__ == '__main__':
  #print clean('*A$C>G<T+23ACACACACACACACACACAGGGGT^JC-2TTC-11AAAAAAAAAAAC+4TTTTC')
  #callReads(*sys.argv[1:])
  binRegions("table.txt")
