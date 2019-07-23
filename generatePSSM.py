from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#Read trimmed alignment
align = AlignIO.read('dF106TrimmedAlignment.fasta', 'fasta')

#UniProtKB/TrEMBL PROTEIN DATABASE RELEASE 2019_05 STATISTICS
aaPerc = {'A':9.17,
          'Q':3.78,
          'L':9.90,
          'S':6.65,
          'R':5.75,
          'E':6.18,
          'K':4.94,
          'T':5.54,
          'N':3.83,
          'G':7.33,
          'M':2.38,
          'W':1.30,
          'D':5.47,
          'H':2.19,
          'F':3.92,
          'Y':2.91,
          'C':1.20,
          'I':5.67,
          'P':4.86,
          'V':6.90}

#Dictionary for the counts of amino acids in the alignement positions
counts = {}

#Loop over all positions in the alignment
for i in range(len(align[0])):
    counts[i] = {}
    
    #Loop over all sequences in the position
    for aa in align[:,i]:
        
        #Count the occurence of each amino acids
        try:
            counts[i][aa] += 1
        except:
            counts[i][aa] = 1

#Convert to a pandas dataframe
df = pd.DataFrame(counts)
df = df.fillna(0)

#Calculate the frequencies with pseudo-counts
PSSM = (df + 1) / (df.sum() + 1)
PSSM = PSSM.T

#Drop information on gaps
PSSM = PSSM.drop(['-'],axis=1)

#Filter out the position with too many gaps
toDrop = []
for i in range(len(align[0])):
    try:
        if counts[i]['-']/len(align) > .9: #Position 85 & 86 out
            toDrop.append(i)
    except:
        pass
print(toDrop)
#No positions in the mutated region are removed
PSSM = PSSM.drop(toDrop) #Remove position with too many gaps

#Calculate the information content by dividing with the uniprot frequencies
PSSMinf = pd.DataFrame()
for aa in aaPerc.keys():
    PSSMinf[aa] = np.log2(PSSM[aa]/(aaPerc[aa]/100))

#Melt the dataframe to be merged with the scores
forPlot = PSSMinf.reset_index()

#Save the PSSM
PSSMinf.to_csv('PSSM.csv')