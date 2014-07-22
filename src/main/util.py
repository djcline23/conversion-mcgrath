import scipy.stats

def cutoffs(binSize):
  """Figure out cutoffs for the given bin size.
  Returns a tuple of cutoffs. Less than the first number
  means the bin is most likely M, greater than the second
  number means the bin is most likely N, in between is 
  heterozygous"""
  
  binSize = int(binSize)  
  #Binomial for number of Ns in bin
  mProb = scipy.stats.binom.pmf(range(binSize+1),binSize,0.01)
  hProb = scipy.stats.binom.pmf(range(binSize+1),binSize,0.5)
  nProb = scipy.stats.binom.pmf(range(binSize+1),binSize,0.99)
  
  mThresh = -1
  nThresh = binSize+1

  for i in range(binSize+1):
    if(hProb[i] > mProb[i]):
      mThresh = i
      break

  for i in range(binSize, -1, -1):
    if(hProb[i] > nProb[i]):
      nThresh = i
      break;

  return float(mThresh), float(nThresh)

if __name__ == '__main__':
  print cutoffs(2)