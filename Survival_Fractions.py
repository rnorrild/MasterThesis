import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import Seq

def color_map(v):
    """ Converts a score (v) into the corresponding color of the heat map.
    The colors are a bit different and the scaling is also different
    """
    bottom = (249/255,59/255,6/255,1)
    top = (6/255,249/255,35/255,0)
    
    if v > 0.5:
        return (38/255,217/255,59/255,float((v-.5)*5 + 0.05))
    else:
        return (249/255,59/255,6/255,float((.5-v)*5 + 0.05))

samples_to_plot = [4,6,7,9]

#Load the dataframes from csv files containing the names of
#substitutions in the different quantile slices created by the 
#script TopScores.py. Another column in the csv specifies the 
#average score of these substitutions
toplot = {}
for sample in samples_to_plot:
    toplot[sample] = pd.read_csv(str(sample) + '_ToPlot.csv')


#Load the data and filter for the analysis
sample_data = {}
for i in np.arange(1,10):
    sample_data[i] = pd.read_csv(''.join(['../190612_NewCutOff/Sample',
                                          str(i),
                                          '_Filtered_AA_Counts.csv']))
    
    #Use only sequences in the input lib
    if i >3:
        input_lib_nr = (i-1)% 3 +1
        sample_data[i] = sample_data[i].merge(sample_data[input_lib_nr]
                                              ,on=['Sequences','AA_Sequences']
                                              ,how='inner')    
    
    #WT sequences are saved to seperate dataframe
    wt_sequences = (sample_data[i]
                    .loc[sample_data[i]['AA_Sequences'].isnull()])
    wt_sequences = (wt_sequences
                    .loc[wt_sequences['Sequences']
                         .notnull()]) #No information is in the unmutated sequence
    
    #Manage non-wildtype sequences alone
    sample_data[i] = (sample_data[i]
                      .loc[sample_data[i]['AA_Sequences'].notnull()])
    
    #Only unique amino acid sequences are of interest

    sample_data[i] = (sample_data[i]
                      .drop_duplicates(subset='AA_Sequences'))
    
    #Sequences with stop codons are omitted from this analysis
    sample_data[i] = (sample_data[i]
                      .loc[~(sample_data[i]
                             ['AA_Sequences']
                             .str.contains('\*',na=False))])

        
    
    #WT sequences are re-introduced
    sample_data[i] = sample_data[i].merge(wt_sequences
                                          ,how='outer')

#A loop to add a column in the dataframes with the amount of substitutions
for sample_nr, d in sample_data.items(): 
    #Assing all to 0 so wildtype will be zero
    for index, row in (d['AA_Sequences']
                       .loc[d['AA_Sequences']
                            .isnull()].iteritems()):
        d.at[index,'#ofSubs'] = 0
        
    #Calculate the amount of substitutions for non wildtype
    for index, row in (d['AA_Sequences']
                       .loc[d['AA_Sequences']
                            .notnull()].iteritems()):
        d.at[index,'#ofSubs'] = len(row.split())
        
#Dictionaries to hold the date for each sample
survival_fraction = {}
errors = {}

AddSub_survival_fraction = {}
AddSub_errors = {}

#Calculate the survival fraction
for sample_nr in samples_to_plot:
    input_lib_nr = (sample_nr-1)% 3 +1
    
    #The survival fraction at zero substitutions is set to one with no error
    survival_fraction[sample_nr] = []
    errors[sample_nr] = []
    
    #The plot will go from 0 to 15 substitutions. This loop is to calculate all
    #sequences and the subsection of sequences are calculated in another loop
    #below.
    for muts in np.arange(0,15):
        #Slice the dataset to sequences with the given amount of substitutions
        df = (sample_data[sample_nr].loc[sample_data[sample_nr]['#ofSubs'] == muts])
        df_input = (sample_data[input_lib_nr]
                    .loc[sample_data[input_lib_nr]
                         ['#ofSubs'] == muts])
        
        #Check if there is enogh data to make the plots
        if ((len(df) + len(df_input)) > 15) and (len(df) > 0) :
            
            #A try/except to check if there are any sequences
            try:
                survival_fraction[sample_nr].append(len(df)/len(df_input))
                if len(df)/len(df_input) == 50.0:
                    print(len(df),len(df_input))
            except:
                survival_fraction[sample_nr].append(np.nan)
            
            #Calculate the errors
            #If there are no growing sequences the error is arbitrarily set to 2
            if len(df) == 0:
                errors[sample_nr].append(np.nan)
            else:
                #The fractional standard errors (poisson) are added in quadrature
                #and multiplied by 1.96 to get the 95% CI. This is in principle
                #only true for large numbers when the poisson distribution looks
                #like a normal distribution.
                try:
                    (errors[sample_nr]
                     .append(np.sqrt((np.sqrt(len(df))/len(df))**2
                                     +(np.sqrt(len(df_input))/len(df_input))**2
                                     )* len(df)/len(df_input) * 1.96)) 
                except:
                    errors[sample_nr].append(2)
        
        else:
            survival_fraction[sample_nr].append(np.nan)
            errors[sample_nr].append(np.nan)                
    
    #The subsection graphs are calculated in the following loop    
    AddSub_survival_fraction[sample_nr] = {}
    AddSub_errors[sample_nr] = {}    
    
    #Loop over all the quantile slices loaded in the beginning
    for i in range(len(toplot[sample])):
        AddSub_survival_fraction[sample_nr][i] = []
        AddSub_errors[sample_nr][i] = []         
        
        #The substitutions in the format "L11P|D83V" to use in regex
        add_sub = toplot[sample_nr]['Seqs'][i]
        
        #This is also calculated from 0 to 15 substitutions
        for muts in np.arange(0,15):
            #Slice the dataset to sequences with the given amount of substitutions
            df = (sample_data[sample_nr].loc[sample_data[sample_nr]['#ofSubs'] == muts])
            df_input = (sample_data[input_lib_nr]
                        .loc[sample_data[input_lib_nr]
                             ['#ofSubs'] == muts])
            
            #Check if there is enogh data to make the plots            
            if ((len(df) + len(df_input)) > 25) and (len(df) > 0) :
                
                if not add_sub == '':    
                    #Only use the part of the dataset with the given mutation
                    AddSub_df = (df.loc[df['AA_Sequences']
                                        .str.contains(add_sub,na=False)])
                    AddSub_df_input = (df_input.loc[df_input['AA_Sequences']
                                                    .str.contains(add_sub,na=False)])
                    
                    #A try/except to check if there are any sequences
                    try:
                        AddSub_survival_fraction[sample_nr][i].append(len(AddSub_df) /
                                                                      len(AddSub_df_input))
                    except:
                        AddSub_survival_fraction[sample_nr][i].append(np.nan)
                        
                    #Calculate the errors               
                    #If there are no growing sequences the error is arbitrarily set to 2
                    if len(AddSub_df) == 0:
                        AddSub_errors[sample_nr][i].append(np.nan) 
                    else:
                        #The fractional standard errors (poisson) are added in quadrature
                        #and multiplied by 1.96 to get the 95% CI. This is in principle
                        #only true for large numbers when the poisson distribution looks
                        #like a normal distribution.
                        try:
                            (AddSub_errors[sample_nr][i]
                             .append(np.sqrt((np.sqrt(len(AddSub_df))/len(AddSub_df))**2
                                             +
                                             (np.sqrt(len(AddSub_df_input))/len(AddSub_df_input))**2
                                             ) * len(AddSub_df)/len(AddSub_df_input) * 1.96)) 
                        except:
                            AddSub_errors[sample_nr][i].append(2)
            else:
                
                if not add_sub == '':
                    AddSub_survival_fraction[sample_nr][i].append(np.nan)
                    AddSub_errors[sample_nr][i].append(np.nan)               

########################################################################################
### Plotting 
########################################################################################

title_dict = {4:'[48:72] 30deg',
              6:'[48:97] 30deg',
              7:'[48:72] 37deg',
              9:'[48:97] 37deg',
              1:'[48:72] input library',
              3:'[48:97] input library'}

#Initialize the plot
f, axarr = plt.subplots(2, 2, sharey=False, sharex=False)
f.set_size_inches(10, 10)

for sample_nr in samples_to_plot:
    #Plotting the date in the right subplot
    plot_map = {4:(0,0),
                6:(1,0),
                7:(0,1),
                9:(1,1)}
    subplot_x_coord,subplot_y_coord = plot_map[sample_nr]
    ax = axarr[subplot_y_coord,subplot_x_coord]
    
    #Plot the data for all sequences with error bars
    ax.errorbar(np.arange(0,15),
                survival_fraction[sample_nr],
                errors[sample_nr],
                capsize=5,
                alpha=0.5,
                label = 'All',
                elinewidth = 1)
    
    #Plot the add_sub data
    if not add_sub == '':
        for i in range(len(toplot[sample_nr])):
            ax.plot(np.arange(0,14),
                    AddSub_survival_fraction[sample_nr][i][1:],
                    color=color_map(toplot[sample_nr]['Avgs'][i]),
                    label = round(toplot[sample_nr]['Avgs'][i],3))

    ax.axhline(y=0, color='grey'
               , linestyle='--',alpha=.3)
    ax.axhline(y=1, color='grey'
               , linestyle='--',alpha=.3)   
    
    ax.legend(loc='upper right')
    
    ax.set_xlim([-1,15])    
    ax.set_ylim([-0.15,1.15])
    ax.set_title(title_dict[sample_nr])
    
    ax.set_xlabel('Number of substitutions')
    ax.set_ylabel('Survival fraction')
    
f.tight_layout()
f.show()