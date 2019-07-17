import pandas as pd
import os
import numpy as np
import re
from Bio import Seq
import matplotlib.pyplot as plt

#The cutoffs to check
interval = np.arange(0,100,1)

#The postitions that should not be mutated
notMutated = list(np.arange(61,71)) + [142,143,144,145,146] + list(np.arange(217,227))
notMutated = [str(item) for item in notMutated]

title_dict = {1:'[48:72] input library',
              2:'[74:97] input library',
              3:'[48:97] input library',
              4:'[48:72] 30deg',
              5:'[74:97] 30deg',
              6:'[48:97] 30deg',
              7:'[48:72] 37deg',
              8:'[74:97] 37deg',
              9:'[48:97] 37deg'}

cutoffs = {1:2,
           2:7,
           3:2,
           4:3,
           5:23,
           6:10,
           7:3,
           8:28,
           9:13}

#Load the data that is named like: 'Sample1_AlignmentOutput_DNA_count.csv'
count_files = {}
for file in os.listdir('../190612_NewFormat'):
    if '.csv' in file:
        count_files[file] = pd.read_csv(file)

#Loop through all samples to do the calculations
unMut = {} 
unique_seqs = {}    
for file_name, df in count_files.items():
    unique_seqs[file_name] = []
    unMut[file_name] = [0 for i in interval]
    
    counter = 0
    for cutoff_over in interval:
        #Check how many sequences that are left
        unique_seqs[file_name].append(len(df['Sequences']
                                          .loc[(df['Sequences'].notnull()) 
                                               & (df['Counts'] > cutoff_over)]))
        
        #Loop over the mutations in all sequences that are left.
        #This is could be done more effeciently 
        for index, row in (df['Sequences'].loc[(df['Sequences'].notnull()) & (df['Counts'] > cutoff_over)].iteritems()):
            for sub in row.split():
                if sub[1:-1] in notMutated:
                    unMut[file_name][counter] += 1 
        
        counter += 1    


#Plot all 
f, axarr = plt.subplots(3, 3, sharey=False, sharex=False)
f.set_size_inches(18.5, 10.5)

to_plot = {}   

count = 1
for i in range(len(axarr)):
    for j in range(len(axarr[i])):
        #Get the data from the righ sample
        temp_name = "".join(['Sample',str(count),'_AlignmentOutput_DNA_count.csv'])
        
        #Plot the data as a fraction (illegalMuts/sequences)
        axarr[i][j].scatter(np.log10(np.array(unique_seqs[temp_name]))
                            ,np.array(unMut[temp_name])/np.array(unique_seqs[temp_name])
                            ,marker='.',label='Un-mutated',color='k')
        
        #Plot the cutoffs 
        axarr[i][j].axvline(x=np.log10(unique_seqs[temp_name][cutoffs[count]])
                            ,color='tab:green'
                            ,linestyle = '--'
                            ,alpha = .3)        
        
        #Define the axis
        axarr[i][j].set_ylim(bottom=0)
        axarr[i][j].set_ylim([0,0.15])
        axarr[i][j].set_xlim([1.5,5])

        #Dummy names to make room
        if count == 7:
            axarr[i][j].set_xlabel(' ')
            axarr[i][j].set_ylabel(' ')

        count += 1

f.text(0.5, 0.07, 'log10(Sequences)', ha='center',fontsize=11)
f.text(0.08, 0.5, r'$\frac{IllegalMutations}{TotalSequences}$', va='center', rotation='vertical',fontsize=14)  

