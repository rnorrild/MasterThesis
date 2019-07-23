import pandas as pd
pd.set_option('display.max_columns', 500)
import matplotlib.pyplot as plt
import numpy as np
from Bio.Seq import Seq


#Generate wild type dataframe
template = Seq(''.join(['CCGTAAAGATAACGACGAAGAGGCCAAGAAGGTTGAGTATATT',
                        'GTGCGCGAACTGGCGCAGGAATTTGACGGTCTGATCATGGTTT',
                        'TCGAGCTGGACACGAACAAGGCACCGGAGATCGCGAAAAAGTA',
                        'CAATATCACCACCACCCCGACTGTCGCATTTTTCAAAAATGGC',
                        'GAGGTCAAGAGCGTTCTGATTGGCGCGATTCCAAAAGACCAGC',
                        'TGCGTGATGAAATCCTGAAATATCTGGGTCACCATCATCACCA',
                        'TCACGGTACCAAACCGTACCAACGTCAGTTCATCGA']))

align_temp = template[70:217].translate() #Mutated region with an extra base (70) for translation of

wt_dict = {}

for i in range(len(align_temp)):
    wt_dict[tuple((i+48,align_temp[i]))] = 1

#Convert dict into dataframe
wt_df = pd.DataFrame(wt_dict.items(),columns=['unpack','is_wildtype'])
wt_df[['Position', 'Mutation']] = pd.DataFrame(wt_df['unpack'].
                                                   tolist(), 
                                                   index=wt_df.index)

#Rotate the dataframe to allow for plotting
wt_df = wt_df[['Position','Mutation','is_wildtype']].copy()
wt_df.sort_values(by='Position')


#####################################################################################

#Read the data from gremlin
pssm = pd.read_csv('PSSM.csv',sep=',')
pssm = pssm.rename(index=str, columns={"Unnamed: 0": "Position"})



#The following command collapses the columns into 
#indentifier variable in one single column to allow for merging
#with the scores
pssm = pd.melt(pssm,
               id_vars=['Position'],
               value_vars=pssm.columns[1:],
               var_name='AA',
               value_name='PSSM_score')

#Get the PSSM scores of the wild type
wt_df = wt_df.merge(pssm,
                    how='left',
                    left_on=['Position','Mutation'],
                    right_on=['Position','AA'])

#Load data from 190529_NewHeatmap scores
oligo_4 = pd.read_csv('../190612_NewHeatMap/Oligo4_30deg_data.csv')
oligo_45 = pd.read_csv('../190612_NewHeatMap/Oligo4+5_30deg_data.csv')
oligo_4_37 = pd.read_csv('../190612_NewHeatMap/Oligo4_37deg_data.csv')
oligo_45_37 = pd.read_csv('../190612_NewHeatMap/Oligo4+5_37deg_data.csv')

d = {4:oligo_4,
     6:oligo_45,
     7:oligo_4_37,
     9:oligo_45_37}

#Filter the data
for sample_nr, df in d.items():
    d[sample_nr] = (df.loc[(df['Counts'] > 80) & (df['Position'] != 83)]
                    [['Position','Mutation','dgs','Score']]).copy()

#Merge with the PSSM data frame 
for sample_nr, df in d.items():
    d[sample_nr] = df.merge(pssm,
                            how='inner',
                            left_on=['Position','Mutation'],
                            right_on=['Position','AA'])

#Merge with the wildtype scores 
for sample_nr, df in d.items():
    d[sample_nr] = (df.merge(wt_df,
                             how="left",
                             left_on='Position',
                             right_on='Position',
                             suffixes=['','_wt']))

    #Calculate the normalized scores. Not entirely sure on the minus but seems good
    d[sample_nr]['norm_PSSM_score'] = d[sample_nr]['PSSM_score'] - d[sample_nr]['PSSM_score_wt']

##############################################################################

title_dict = {4:'[48:72] 30deg',
              6:'[48:97] 30deg',
              7:'[48:72] 37deg',
              9:'[48:97] 37deg',
              1:'[48:72] input library',
              3:'[48:97] input library'}


#Initialise the figure
f, axarr = plt.subplots(1,4, sharey=False)
f.set_size_inches(13, 3.3)

for i in [0,1,2,3]:
    #Set the axis
    axarr[i].set_xlim([0,1])
    axarr[i].set_ylim([-12,12])
    
    #Plot the lines representing neutral
    axarr[i].axvline(x=0.5
                ,color='gray'
                ,linestyle='--'
                ,alpha = 0.5)
    axarr[i].axhline(y=0
                ,color='gray'
                ,linestyle='--'
                ,alpha = 0.5)
    
    #Annotate
    axarr[i].set_xlabel('Score')
    axarr[i].set_ylabel('PSSM substitution score')

i = 0
for sample_nr, df in d.items():
    #Loop over the four samples and plot them in the same figure
    axarr[i].set_title(title_dict[sample_nr])
    
    #Plotting the [48:72] region
    axarr[i].scatter(df.loc[df['Position'] < 73]['Score']
                     ,df.loc[df['Position'] < 73]['norm_PSSM_score']
                     ,marker='.'
                     ,s = 8
                     ,label='[48:72] region')
    
    #Plotting the [74:97] region
    axarr[i].scatter(df.loc[df['Position'] > 73]['Score']
                     ,df.loc[df['Position'] > 73]['norm_PSSM_score']
                     ,marker='.'
                     ,s = 8
                     #,color='red'
                     ,label='[74:97] region')
    
    #Calculate and annotate the pearson correlation coefficient
    axarr[i].text(.7,-4,r'$\rho$=' + str(round(df[['Score','norm_PSSM_score']]
                                               .corr(method='pearson')['Score']['norm_PSSM_score'],4)))

    if i ==1:
        axarr[i].legend(loc='upper right')
    i += 1
f.tight_layout()

