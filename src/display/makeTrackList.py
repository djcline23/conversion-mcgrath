#!/usr/bin/env python                                                                  

'''                                                                                     
djc 7/21/14                                                                         
'''

import sys
import os

readColors = {"U" : "gray", "H" : "purple", "M" : "red", "N" : "blue"}

def writeStart(outFile):
  """Write out the constant start of trackList.json"""
  outFile.write('{\n')
  outFile.write('   "tracks" : [\n')
  outFile.write('      {\n')
  outFile.write('         "chunkSize" : 20000,\n')
  outFile.write('         "storeClass" : "JBrowse/Store/Sequence/StaticChunked",\n')
  outFile.write('         "urlTemplate" : "seq/{refseq_dirpath}/{refseq}-",\n')
  outFile.write('         "category" : "Reference sequence",\n')
  outFile.write('         "type" : "SequenceTrack",\n')
  outFile.write('         "label" : "DNA",\n')
  outFile.write('         "key" : "Reference sequence"\n')
  outFile.write('      },\n')
  
def writeEnd(outFile):
  outFile.write('   ],\n   "formatVersion" : 1\n}')

def writePred(pred, outFile):
  pred = pred.replace('.gff', '')
  
  outFile.write('      {\n')
  outFile.write('         "style" : {\n')
  outFile.write('            "className" : "feature"\n')
  outFile.write('         },\n')
  outFile.write('         "key" : "' + pred + '",\n')
  outFile.write('         "storeClass" : "JBrowse/Store/SeqFeature/NCList",\n')
  outFile.write('         "trackType" : "CanvasFeatures",\n')
  outFile.write('         "urlTemplate" : "tracks/' + pred + '/{refseq}/trackData.json",\n')
  outFile.write('         "compress" : 0,\n')
  outFile.write('         "type" : "CanvasFeatures",\n')
  outFile.write('         "label" : "' + pred + '",\n')
  outFile.write('         "style" : {\n')
  outFile.write('            "color" : function(feature) {\n')
  outFile.write('               var pred = feature.get("prediction");\n')
  outFile.write('               if(pred == "M") return "red";\n')
  outFile.write('               else if(pred == "N") return "blue";\n')
  outFile.write('               else if(pred == "H") return "purple";\n')
  outFile.write('               else return "gray";\n')
  outFile.write('            }\n')
  outFile.write('         }\n')
  outFile.write('      },\n')
  
def writeRead(read, outFile):
  read = read.replace('.gff', '')
  
  outFile.write('      {\n')
  outFile.write('         "style" : {\n')
  outFile.write('            "className" : "feature",\n')
  outFile.write('            "color" : "' + readColors[read[-1]] + '"\n')
  outFile.write('         },\n')
  outFile.write('         "key" : "' + read + '",\n')
  outFile.write('         "storeClass" : "JBrowse/Store/SeqFeature/NCList",\n')
  outFile.write('         "trackType" : "CanvasFeatures",\n')
  outFile.write('         "urlTemplate" : "tracks/' + read + '/{refseq}/trackData.json",\n')
  outFile.write('         "compress" : 0,\n')
  outFile.write('         "type" : "CanvasFeatures",\n')
  outFile.write('         "label" : "' + read + '",\n')
  outFile.write('         "displayMode" : "collapsed",\n')
  outFile.write('         "maxFeatureScreenDensity" : 10000\n')
  outFile.write('      },\n')

def makeTrackList():
  """Makes the tracklist file for JBrowse from the reads and preds directories"""
  with open('trackList.json', 'w') as outFile:
    writeStart(outFile)
    for pred in os.listdir('preds'):
      writePred(pred, outFile)
    for read in os.listdir('reads'):
      writeRead(read, outFile)
    writeEnd(outFile)

if __name__ == '__main__':
  makeTrackList(*sys.argv[1:])
