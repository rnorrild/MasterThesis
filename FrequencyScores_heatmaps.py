import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy
from Bio import Seq

#Read in data from csv files
rosetta = pd.read_csv('../190320_RosettaDDGs/Norm_Rosetta_ddg.csv',sep='\t',index_col = 0)
rosetta['Position'] = rosetta['Position'] - 1 #Correction of the position numbering

#Read files containg the DNA and corresponding amino acid sequences and their read counts
sample_data = {}
for i in np.arange(1,10):
    sample_data[i] = pd.read_csv(''.join(['../190612_NewCutOff/Sample',
                                          str(i),
                                          '_Filtered_AA_Counts.csv']))
    
    #No wild type amino acid sequences
    sample_data[i] = sample_data[i].loc[sample_data[i]['AA_Sequences'].notnull()]
    
    #Only unique amino acid sequences 
    sample_data[i] = sample_data[i].drop_duplicates(subset='AA_Sequences')
    
#Samples to plot
samples_to_plot = [4,6,7,9]

#Order of amino acids in the heatmaps
aa_order = ['G','A','V','L','M','I','F','Y','W','K','R',
            'H','D','E','S','T','C','N','Q','P','*','Hydrophobic','Hydrophilic']

#Dictionary to map the input libray
input_lib_map = {4:1,
                 5:2,
                 6:3,
                 7:1,
                 8:2,
                 9:3}

#Perform the calculations for each of the selected datasets

#Dictionaries to use outside of the loop
input_counts_df = {}
growing_counts_df = {}
non_growing_counts_df = {}

len_growing = {}
len_non_growing = {}

for sample_nr in samples_to_plot:
    #Dictionaries to contain the counts before making the dataframes
    growing_count = {}
    non_growing_count = {}
    input_count = {}
    
    #Defining the dataset to use
    growing = sample_data[sample_nr]
    input_lib = sample_data[input_lib_map[sample_nr]]

    #Define the dataset for the nongrowing variants
    non_growing = input_lib.loc[~input_lib['AA_Sequences']
                                .isin(growing['AA_Sequences'])]   
    
    #For calculating the scores later. The number of total sequences
    len_growing[sample_nr] = len(growing)
    len_non_growing[sample_nr] = len(non_growing)  
    
    
    #A loop for counting each occurence of a Substitution in the grwoing dataset
    for index, row in growing['AA_Sequences'].iteritems():
        for mut in row.split(): #Loop over all Substitution in the sequence  
            pos = int(mut[1:-1])
            new = mut[-1]
            try:
                growing_count[tuple((pos,new))] += 1
            except:    
                growing_count[tuple((pos,new))] = 1

    
    #The non_growing sample is counted in the same fashion as the growing sample
    for index, row in non_growing['AA_Sequences'].iteritems():
        for mut in row.split(): #Loop over all Substitution in the sequence  
            pos = int(mut[1:-1])
            new = mut[-1]
            try:
                non_growing_count[tuple((pos,new))] += 1
            except:    
                non_growing_count[tuple((pos,new))] = 1
    
    #A loop for counting the number of occurences of all Substitutions in the input library
    for index, row in input_lib['AA_Sequences'].iteritems():
        for mut in row.split(): #Loop over all Substitution in the sequence  
            pos = int(mut[1:-1])
            new = mut[-1]
            try:
                input_count[tuple((pos,new))] += 1
            except:    
                input_count[tuple((pos,new))] = 1
    
    #These lines makes a dataframe for each of the samples gathered in a dictionary
    #Make dataframe from the input_counts_df
    input_counts_df[sample_nr] = (pd.DataFrame(input_count.items(),
                                               columns=['unpack','Counts']))
    
    (input_counts_df[sample_nr]
     [['Position', 'Substitution']]) = pd.DataFrame(input_counts_df[sample_nr]
                                                    ['unpack']
                                                    .tolist(), 
                                                    index=input_counts_df[sample_nr].index)


    #Make dataframe from the growing_counts_df
    growing_counts_df[sample_nr] = (pd.DataFrame(growing_count.items(),
                                                 columns=['unpack','Counts']))
    
    (growing_counts_df[sample_nr]
     [['Position', 'Substitution']]) = pd.DataFrame(growing_counts_df[sample_nr]
                                                    ['unpack']
                                                    .tolist(), 
                                                    index=growing_counts_df[sample_nr].index)
        
    #Make dataframe from the non_growing_counts_df
    non_growing_counts_df[sample_nr] = (pd.DataFrame(non_growing_count.items(),
                                                     columns=['unpack','Counts']))
    
    (non_growing_counts_df[sample_nr]
     [['Position', 'Substitution']]) = pd.DataFrame(non_growing_counts_df[sample_nr]
                                                    ['unpack']
                                                    .tolist(), 
                                                    index=non_growing_counts_df[sample_nr].index)


final_dfs = {}
#Merge the counts data to the Rosetta dataframes
for sample_nr in samples_to_plot:
    final_dfs[sample_nr] = (input_counts_df[sample_nr]
                            .merge(rosetta,
                                   how='left', 
                                   left_on=['Position','Substitution'], 
                                   right_on = ['Position','Mutation']))
    
    final_dfs[sample_nr] = pd.merge(final_dfs[sample_nr],
                                    growing_counts_df[sample_nr],
                                    how='left', 
                                    left_on=['Position','Substitution'], 
                                    right_on = ['Position','Substitution'],
                                    suffixes=('','_growing'))    
    
    final_dfs[sample_nr] = pd.merge(final_dfs[sample_nr],
                                    non_growing_counts_df[sample_nr],
                                    how='left', 
                                    left_on=['Position','Substitution'], 
                                    right_on = ['Position','Substitution'],
                                    suffixes=('','_non_growing'))
    
#Scoring the counts
for sample_nr in samples_to_plot:
    df = final_dfs[sample_nr].copy(deep=True)
   
    df['Score'] = (df['Counts_growing'].divide(len_growing[sample_nr])
                   / (df['Counts_growing'].divide(len_growing[sample_nr])
                      + df['Counts_non_growing'].divide(len_non_growing[sample_nr]))) 
    final_dfs[sample_nr] = df



##################
#Generate wild type dataframe
template = Seq.Seq(''.join(['CCGTAAAGATAACGACGAAGAGGCCAAGAAGGTTGAGTATATT',
                            'GTGCGCGAACTGGCGCAGGAATTTGACGGTCTGATCATGGTTT',
                            'TCGAGCTGGACACGAACAAGGCACCGGAGATCGCGAAAAAGTA',
                            'CAATATCACCACCACCCCGACTGTCGCATTTTTCAAAAATGGC',
                            'GAGGTCAAGAGCGTTCTGATTGGCGCGATTCCAAAAGACCAGC',
                            'TGCGTGATGAAATCCTGAAATATCTGGGTCACCATCATCACCA',
                            'TCACGGTACCAAACCGTACCAACGTCAGTTCATCGA']))
align_temp = template[70:217].translate() #Mutated region with an extra base (70) for translation

wt_dict = {}

for i in range(len(align_temp)):
    wt_dict[tuple((i+48,align_temp[i]))] = 1

#Convert dict into dataframe
wt_df = pd.DataFrame(wt_dict.items(),columns=['unpack','is_wildtype'])
wt_df[['Position', 'Substitution']] = pd.DataFrame(wt_df['unpack'].
                                                   tolist(), 
                                                   index=wt_df.index)

#Rotate the dataframe to allow for plotting
wt_df = wt_df[['Position','Substitution','is_wildtype']].copy()
wt_df.sort_values(by='Position')

#Export the WT dataframe
wt_df.to_csv('WT_df.csv')

#Pivot for plotting
wt_df = wt_df.pivot(index='Position', columns = 'Substitution')
wt_df.columns = wt_df.columns.droplevel(0)

wt_df = wt_df.T

##################



################## Plotting all the data as heatmaps
#Sort the wildtype on the aa order
wt_df = wt_df.reindex(aa_order)

#Initialize the plot
f, axarr = plt.subplots(4, 1, sharey=False, sharex=False)
f.set_size_inches(7, 10)

title_dict = {4:'[48:72] 30deg',
              5:'Oligo 5 30deg',
              6:'[48:97] 30deg',
              7:'[48:72] 37deg',
              8:'Oligo 5 37deg',
              9:'[48:97] 37deg'}

count = 0
for sample_nr in [4,7,6,9]:
    subplot_y_coord = count

    #Copy and otate the dataframe to allow for a heatmap
    df1 = (final_dfs[sample_nr]
           .loc[final_dfs[sample_nr]['Counts'] > 80]
           [['Position','Substitution','Score']].copy(deep=True))  
    
    #Filling the empty space in the plot that shifts everything
    for i in np.arange(48,97):
        if not i in list(df1['Position']):
            df1 = df1.append(pd.DataFrame([[i,'*',np.nan]],
                                          columns=['Position','Substitution','Score']))
            
    #Sort the data
    df1.sort_values(by='Position')
    
    #Calculating averages for the plot
    hydrophobic = ['G','A','V','L','M','I','P','F','W']           
    hydrophilic = ['Y','K','R','H','D','E','S','T','C','N','Q',]
    
    #Only the rows where the substitutions are hydrophilic or ydrophilic
    philic = df1.loc[df1['Substitution'].isin(hydrophilic)]
    phobic = df1.loc[df1['Substitution'].isin(hydrophobic)]
    
    #Calculate the average score for each amino acid
    philic = philic.sort_values('Position').groupby('Position')['Score'].mean().reset_index()
    phobic = phobic.sort_values('Position').groupby('Position')['Score'].mean().reset_index()
    
    #Name the columns
    philic['Substitution'] = 'Hydrophilic'
    phobic['Substitution'] = 'Hydrophobic'
    
    #Append the average values to the other data
    df1 = df1.append(philic, sort=False)
    df1 = df1.append(phobic, sort=False)
    
    #Pivot the dataframe for plotting 
    df1 = df1.pivot(index='Position', columns = 'Substitution')
    df1.columns = df1.columns.droplevel(0)
    df1 = df1.T
    
    #Sort the dataframe on the amino acid order
    df1 = df1.reindex(aa_order)
    
    #Round off to two decimals 
    df1 = df1.round(2)
    
    #Plotting the counts
    sns.heatmap(df1,
                cmap=sns.diverging_palette(13, 127, l=50, s=95, n = 1000),
                center = 0.5,
                ax=axarr[subplot_y_coord],
                xticklabels=True, yticklabels=True,
                annot=True,
                annot_kws={"size": 3})
    
    axarr[subplot_y_coord].tick_params(axis='both', which='major', labelsize=5)
    
    #Plotting wild type
    sns.heatmap(wt_df,
                cmap=sns.diverging_palette(127, 10, l=0, s=95, n = 1000),
                cbar=False,
                ax=axarr[subplot_y_coord])

    axarr[subplot_y_coord].set_title(title_dict[sample_nr])
    f.tight_layout()
    
    count += 1


"""
#All correlations plotted (scatter)
cutOff = 80
toPlot = [(4,7),(4,6),(4,9)
          ,(6,9),(7,9),(6,7)]

f, axarr = plt.subplots(2, 3, sharey=False, sharex=False)
f.set_size_inches(18.5, 10.5)

def traverse(item):
    try:
        for i in iter(item):
            for j in traverse(i):
                yield j
    except TypeError:
        yield item

count = 0
for ax in traverse(axarr):
    sample1, sample2 = toPlot[count]
    #All Substitutions
    test = pd.merge(final_dfs[sample1]
                    ,final_dfs[sample2]
                    ,how='inner' 
                    ,on=['Position','Substitution'])
    ax.scatter(test['Score_x'],test['Score_y']
               ,color='gray',alpha=.1,marker='.') 
    
    #Only above the cutoff
    test = pd.merge(final_dfs[sample1].loc[final_dfs[sample1]['Counts'] > cutOff]
                    ,final_dfs[sample2].loc[final_dfs[sample2]['Counts'] > cutOff]
                    ,how='inner' 
                    ,on=['Position','Substitution'])
    ax.scatter(test['Score_x'],test['Score_y'],alpha=.5,color='tab:blue',marker='.')    
    
    correl = (test[['Score_x','Score_y']].corr()['Score_x']['Score_y'])
    ax.text(.7,.2,r'$\rho$=' + str(round(correl,3)))
    
    ax.set_ylim([0,1])
    ax.set_xlim([0,1])
    
    ax.set_xlabel(title_dict[sample1])
    ax.set_ylabel(title_dict[sample2])
    count += 1
"""

"""
#Rosetta correlations for all datasets
for sample in [4,6,7,9]:
        
    
    plt.scatter(final_dfs[sample].loc[(final_dfs[sample]['Counts']>80) & ~(final_dfs[sample]['Position'] == 83)]['Score']
                ,final_dfs[sample].loc[(final_dfs[sample]['Counts']>80) & ~(final_dfs[sample]['Position'] == 83)]['dgs']
                ,marker='.')
    plt.axvline(x=0.5
                ,color='gray'
                ,linestyle='--'
                ,alpha = 0.5)
    plt.axhline(y=0
                ,color='gray'
                ,linestyle='--'
                ,alpha = 0.5)
    plt.ylabel('Rosetta ddg')
    plt.xlabel(title_dict[sample])
    #plt.title('Rosetta cartesian ddg correlation')
    plt.ylim([-12,27])
    plt.xlim([0,.7])
    
    
    correl = (final_dfs[sample].loc[(final_dfs[sample]['Counts']>80) & ~(final_dfs[sample]['Position'] == 83)]
              [['Score','dgs']].corr(method='spearman')['Score']['dgs'])
    plt.text(.1,-8,r'$\rho$=' + str(round(correl,3)))
    
    plt.savefig('RosettaSpearmanCorrelation_Sample' + str(sample) + '.png')
    plt.clf()
"""


"""
#Cutoff plot showing the effect on the correlation of different cutoffs 
cutOff = 80
toPlot = [(4,7),(4,6),(4,9)
          ,(6,9),(7,9),(6,7)]


f, axarr = plt.subplots(2, 3, sharey=True, sharex=False)
f.set_size_inches(18.5, 10.5)

def traverse(item):
    try:
        for i in iter(item):
            for j in traverse(i):
                yield j
    except TypeError:
        yield item

count = 0
for ax in traverse(axarr):
    sample1, sample2 = toPlot[count]
    
    corr_Sp = []
    corr_P = []
    len_merge = []
    for i in np.arange(1,400):
        test = pd.merge(final_dfs[sample1].loc[final_dfs[sample1]['Counts'] > i]
                        ,final_dfs[sample2].loc[final_dfs[sample2]['Counts'] > i]
                        ,how='inner' 
                        ,on=['Position','Substitution'])                    
        
        x = np.array(test['Score_x'])
        y = np.array(test['Score_y'])
        
        new_x = []
        new_y = []
        for i in range(len(x)):
            if not (np.isnan(x[i]) or np.isnan(y[i])):
                new_x.append(x[i])
                new_y.append(y[i])    
            
        c, p = scipy.stats.mstats.spearmanr(new_x,new_y)
        corr_Sp.append(c)
        corr_P.append(np.corrcoef(new_x,new_y)[0,1])
        len_merge.append(len(new_x))
    
    ax.plot(corr_P, color='k', label='Pearsons')
    ax.plot(corr_Sp, color='grey', label='Spearman')
    
    ax.legend()
    ax2 = ax.twinx()
    ax2.plot(len_merge, color='tab:blue')
    ax2.tick_params('y',labelcolor='tab:blue')        
    
    ax.set_title(''.join([title_dict[sample1],' vs ',title_dict[sample2]]))
    count += 1

f.text(0.5, 0.06, 'Cutoff', ha='center',fontsize=11)
f.text(0.08, 0.5, 'Spearman correlation coefficient'
       , va='center', rotation='vertical',fontsize=11)  
f.text(0.94, 0.5, 'Number of substitutions'
       , va='center', rotation='vertical',fontsize=11,color='tab:blue')
f.suptitle('Cutoff effect on correlation and number of substitutions')
"""
"""
#All correlations plotted against the number of substitutions
cutOff = 80
toPlot = [(4,7),(4,6),(4,9)
          ,(6,9),(7,9),(6,7)]


f, axarr = plt.subplots(2, 3, sharey=True, sharex=False)
f.set_size_inches(18.5, 10.5)

def traverse(item):
    try:
        for i in iter(item):
            for j in traverse(i):
                yield j
    except TypeError:
        yield item

count = 0
for ax in traverse(axarr):
    sample1, sample2 = toPlot[count]
    
    corr = []
    len_merge = []
    for i in np.arange(1,600):
        test = pd.merge(final_dfs[sample1].loc[final_dfs[sample1]['Counts'] > i]
                        ,final_dfs[sample2].loc[final_dfs[sample2]['Counts'] > i]
                        ,how='inner' 
                        ,on=['Position','Substitution'])                    
        
        x = np.array(test['Score_x'])
        y = np.array(test['Score_y'])
        
        new_x = []
        new_y = []
        for i in range(len(x)):
            if not (np.isnan(x[i]) or np.isnan(y[i])):
                new_x.append(x[i])
                new_y.append(y[i])    
        
        c, p = scipy.stats.mstats.spearmanr(new_x,new_y)
        corr.append(c)
        #corr.append(np.corrcoef(new_x,new_y)[0,1])
        len_merge.append(len(new_x))
    
    ax.scatter(len_merge,corr, color='tab:blue', marker = '.', alpha = 0.3)
           
    
    ax.set_title(''.join([title_dict[sample1],' vs ',title_dict[sample2]]))
    count += 1

f.text(0.5, 0.06, 'Number of substitutions', ha='center',fontsize=11)
f.text(0.08, 0.5, 'Spearman correlation coefficient'
       , va='center', rotation='vertical',fontsize=11)  

f.suptitle('Correlations as a function of number of substitutions included')
"""