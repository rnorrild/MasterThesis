import pandas as pd
import os
import numpy as np
from Bio import Seq
from matplotlib_venn import venn3
import matplotlib.pyplot as plt

template = Seq.Seq(''.join(['CCGTAAAGATAACGACGAAGAGGCCAAGAAGGTTGAGTATATT',
                            'GTGCGCGAACTGGCGCAGGAATTTGACGGTCTGATCATGGTTT',
                            'TCGAGCTGGACACGAACAAGGCACCGGAGATCGCGAAAAAGTA',
                            'CAATATCACCACCACCCCGACTGTCGCATTTTTCAAAAATGGC',
                            'GAGGTCAAGAGCGTTCTGATTGGCGCGATTCCAAAAGACCAGC',
                            'TGCGTGATGAAATCCTGAAATATCTGGGTCACCATCATCACCA',
                            'TCACGGTACCAAACCGTACCAACGTCAGTTCATCGA']))
align_temp = template[70:217].translate() #Mutated region with an extra base (70) for translation of the first codon 

def generate_mutString(s):
    """ The mutation has to be put back into the full sequence to allow
    for correct translation of the codons
    """
    test_seq = Seq.MutableSeq(str(template))
    #Introduce the mutation in a test string
    for mut in s.split():
        pos = int(mut[1:-1]) - 1 #Numbering from 0 in strings 
        old = mut[0]
        new = mut[-1]
        if old == test_seq[pos]:
            test_seq[pos] = new 
        else:
            print('Initial mutation didnt match')
    
    return test_seq

def get_AA_subs(s):
    """ The MutableSeq object is compared to the align_temp sequence
    to identify substitutions
    """
    test_seq = s.toseq()[70:217].translate() #Translate the mutated region
    substitutions = []
    
    for i in range(len(test_seq)):
        if test_seq[i] != align_temp[i]:
            substitutions.append(''.join([str(align_temp[i]),
                                          str(i+48),
                                          str(test_seq[i]),
                                          ' ']))
    
    return ''.join(substitutions).strip()

def do(s):
    """ Workaround to call two nested functions inside a pandas dataframe
    """
    return get_AA_subs(generate_mutString(s)) 

if __name__ == "__main__":
    #Load the different counts as pandas dataframes
    count_files = {}
    for file in os.listdir('../190612_NewFormat'):
        if '.csv' in file:
            count_files[file] = pd.read_csv(''.join(['../190612_NewFormat/',file]),
                                            skiprows=1,
                                            header=None,
                                            names = ["Sequences", "Counts"])
    
    #Cut-offs chosen from the cutoff plot
    cut_offs = {1:2,
                2:7,
                3:2,
                4:3,
                5:23,
                6:10,
                7:3,
                8:28,
                9:13}
    
    #Filtering the reads based on the cutoffs 
    filtered = {}
    for i in range(1,10):
        name = 'Sample' + str(i) + '_AlignmentOutput_DNA_count.csv'
        filtered[name]  = count_files[name].loc[count_files[name]['Counts'] > cut_offs[i]]
    

    #Loop for translating the mutations and saving the the dataframes
    for file_name, df in filtered.items():
        df['AA_Sequences'] = (df['Sequences'].
                              loc[df['Sequences'].notnull()].
                              apply(do))
        df.to_csv(''.join([file_name[0:8],'Filtered_AA_Counts.csv']))
    
    
    #Plotting the sequences as a Venn diagram and saving it
    title_dict = {1:'[48:72]',
                  2:'[74:97]',
                  3:'[48:97]'}
     
    #Make sets from all the sequences in the different samples
    sets = {}
    for i in np.arange(1,10):
        sets[i] = set((filtered[''.join(['Sample',
                                     str(i),
                                     '_AlignmentOutput_DNA_count.csv'])]['Sequences']))
    
    #Loop over the three input libraries
    for i in np.arange(1,4):             
        venn3([len(sets[i] - sets[i+3] - sets[i+6]),#Input, not in either
               len(sets[i+3] - sets[i] - sets[i+6]),
               len(sets[i].intersection(sets[i+3])- sets[i+6]),
               len(sets[i+6] - sets[i] - sets[i+3]),
               len(sets[i+6].intersection(sets[i]) - sets[i+3]),
               len(sets[i+6].intersection(sets[i+3]) - sets[i]),
               len(sets[i+6].intersection(sets[i],sets[i+3]))],
              set_labels = (['Input library', '30 deg', '37 deg']),
              set_colors = ['g','tab:orange','r'])
        plt.title([title_dict[i]])
        
        
        #Calcul
        plt.text(-1.,-0.55
                 ,''.join(['30 deg overlap fraction: '
                           ,str(round(len(sets[i]
                                          .intersection(sets[i+3]))
                                      /len(sets[i+3]),3))]))
        plt.text(-1.,-0.6
                 ,''.join(['37 deg overlap fraction: '
                           ,str(round(len(sets[i]
                                          .intersection(sets[i+6]))
                                      /len(sets[i+6]),3))]))
        plt.text(-1.,-0.65
                 ,''.join(['Survival fraction at 30 deg: '
                           ,str(round(len(sets[i+3])
                                      /len(sets[i]),3))]))
        plt.text(-1.,-0.7
                 ,''.join(['Survival fraction at 37 deg: '
                           ,str(round(len(sets[i+6])
                                      /len(sets[i]),3))]))
        plt.text(-1.,-0.75
                 ,''.join(['37 deg / 30 deg overlap fraction: '
                           ,str(round(len(sets[i+3]
                                          .intersection(sets[i+6]))
                                      /len(sets[i+6]),3))]))        
        plt.savefig(''.join([str(i),'_Venn','.png'])) 
        plt.clf()
        