# this notebook adds the sex, jackpot status to each sample 

import networkx as nx
import numpy as np
import matplotlib.pylab as plt
import random
from scipy.stats import ttest_ind
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
from networkx.drawing.nx_agraph import graphviz_layout

def normalize_key(key):
    return tuple(sorted(key))



#load READ result files 
file_1= '/gpfs/scratch/akwakye/AIM_1/lcmlkin/SC_13_14_READ_results'
file_2='/gpfs/scratch/akwakye/AIM_1/lcmlkin/SC_13_14_norm_meansP0_AncientDNA_normalized'
READ_files=open(file_1,'r')
READ_list=READ_files.read()
SC_13_14_READ_result=READ_list.split('\n')
if SC_13_14_READ_result[-1] =='':
    del(SC_13_14_READ_result[-1])

READ_files=open(file_2,'r')
READ_list=READ_files.read()
SC_13_14_READ_values=READ_list.split('\n')
if SC_13_14_READ_values[-1] =='':
    del(SC_13_14_READ_values[-1])

    
    
file_1= '/gpfs/scratch/akwakye/AIM_1/lcmlkin/SC13_vs_SC15_READ_results'
file_2='/gpfs/scratch/akwakye/AIM_1/lcmlkin/SC13_vs_SC15_meansP0_AncientDNA_normalized'
READ_files=open(file_1,'r')
READ_list=READ_files.read()
SC13_vs_SC15_READ_result=READ_list.split('\n')
if SC13_vs_SC15_READ_result[-1] =='':
    del(SC13_vs_SC15_READ_result[-1])
del(SC13_vs_SC15_READ_result[0])

READ_files=open(file_2,'r')
READ_list=READ_files.read()
SC13_vs_SC15_READ_values=READ_list.split('\n')
if SC13_vs_SC15_READ_values[-1] =='':
    del(SC13_vs_SC15_READ_values[-1])

del(SC13_vs_SC15_READ_values[0])    

    
file_1= '/gpfs/scratch/akwakye/AIM_1/lcmlkin/SC13_vs_SC17_READ_results'
file_2='/gpfs/scratch/akwakye/AIM_1/lcmlkin/SC13_vs_SC17_meansP0_AncientDNA_normalized'
READ_files=open(file_1,'r')
READ_list=READ_files.read()
SC13_vs_SC17_READ_result=READ_list.split('\n')
if SC13_vs_SC17_READ_result[-1] =='':
    del(SC13_vs_SC17_READ_result[-1])
del(SC13_vs_SC17_READ_result[0])

READ_files=open(file_2,'r')
READ_list=READ_files.read()
SC13_vs_SC17_READ_values=READ_list.split('\n')
if SC13_vs_SC17_READ_values[-1] =='':
    del(SC13_vs_SC17_READ_values[-1])

del(SC13_vs_SC17_READ_values[0])    



file_1= '/gpfs/scratch/akwakye/AIM_1/lcmlkin/SC14_vs_SC15_READ_results'
file_2='/gpfs/scratch/akwakye/AIM_1/lcmlkin/SC14_vs_SC15_meansP0_AncientDNA_normalized'
READ_files=open(file_1,'r')
READ_list=READ_files.read()
SC14_vs_SC15_READ_result=READ_list.split('\n')
if SC14_vs_SC15_READ_result[-1] =='':
    del(SC14_vs_SC15_READ_result[-1])
del(SC14_vs_SC15_READ_result[0])

READ_files=open(file_2,'r')
READ_list=READ_files.read()
SC14_vs_SC15_READ_values=READ_list.split('\n')
if SC14_vs_SC15_READ_values[-1] =='':
    del(SC14_vs_SC15_READ_values[-1])

del(SC14_vs_SC15_READ_values[0])    




file_1= '/gpfs/scratch/akwakye/AIM_1/lcmlkin/SC14_vs_SC17_READ_results'
file_2='/gpfs/scratch/akwakye/AIM_1/lcmlkin/SC14_vs_SC17_meansP0_AncientDNA_normalized'
READ_files=open(file_1,'r')
READ_list=READ_files.read()
SC14_vs_SC17_READ_result=READ_list.split('\n')
if SC14_vs_SC17_READ_result[-1] =='':
    del(SC14_vs_SC17_READ_result[-1])
del(SC14_vs_SC17_READ_result[0])

READ_files=open(file_2,'r')
READ_list=READ_files.read()
SC14_vs_SC17_READ_values=READ_list.split('\n')
if SC14_vs_SC17_READ_values[-1] =='':
    del(SC14_vs_SC17_READ_values[-1])

del(SC14_vs_SC17_READ_values[0])    


file_1= '/gpfs/scratch/akwakye/AIM_1/lcmlkin/SC15_vs_SC17_READ_results'
file_2='/gpfs/scratch/akwakye/AIM_1/lcmlkin/SC15_vs_SC17_meansP0_AncientDNA_normalized'
READ_files=open(file_1,'r')
READ_list=READ_files.read()
SC15_vs_SC17_READ_result=READ_list.split('\n')
if SC15_vs_SC17_READ_result[-1] =='':
    del(SC15_vs_SC17_READ_result[-1])
del(SC15_vs_SC17_READ_result[0])

READ_files=open(file_2,'r')
READ_list=READ_files.read()
SC15_vs_SC17_READ_values=READ_list.split('\n')
if SC15_vs_SC17_READ_values[-1] =='':
    del(SC15_vs_SC17_READ_values[-1])

del(SC15_vs_SC17_READ_values[0])    


READ_values=SC_13_14_READ_values+SC13_vs_SC15_READ_values+SC13_vs_SC17_READ_values+SC14_vs_SC15_READ_values+SC14_vs_SC17_READ_values+ SC15_vs_SC17_READ_values  

READ_result=SC_13_14_READ_result+SC13_vs_SC15_READ_result+SC13_vs_SC17_READ_result+SC14_vs_SC15_READ_result+SC14_vs_SC17_READ_result+SC15_vs_SC17_READ_result

#load males 
Male_files=open('All_males.list','r')
Male_list=Male_files.read()
Males=Male_list.split('\n')
if Males[-1] =='':
    del(Males[-1])

Males_all_samps=[]
Males_all_samps_dict={}
for i in range(len(Males)):
    k=Males[i].split('/')[-1].split('.')[0].split('_')[0]
    Males_all_samps.append(k)
    Males_all_samps_dict[k]=k

#load jackpot carriers 
file=open('jackpots.txt','r') 
File_list=file.read()
jackpots_list=File_list.split('\n')
if jackpots_list[-1] =='':
    del(jackpots_list[-1])
All_jackpots=jackpots_list

    
Jackpot='Jackpot'
Non_jackpot='Non_jackpot'
READ_result_merged_jack_sex=[]
for i in range(1,len(READ_result)):
    k=READ_result[i].split('\t')
    kk=split_string(k[0])
    #kkk=split_string(k[0])[1]
    year1=kk[0].split('-')[1]
    year2=kk[1].split('-')[1]
    kkk=[year1, year2]
    m=READ_values[i].split(' ')
    #if  k[1] =='Second Degree':
    
    if  k[1] =='First Degree' or k[1] =='Second Degree' or k[1].split('/')[0] =='IdenticalTwins':  
        #if all(e in Males_all_samps for e in kk):          
        if any(e in Males_all_samps for e in kk): #one male and one female
            index_male=[iii for iii, e in enumerate(kk) if e in Males_all_samps]
            if len(index_male) == 1: 
                if index_male[0] == 0: #first is male and second is female 
                    
                    if all(e in All_jackpots for e in kk): # both samples are jackpots 
                        out_f=kk[0], kk[1], k[1], m[3],  m[1], 'M', 'F', Jackpot, Jackpot, kkk[0], kkk[1]
                        READ_result_merged_jack_sex.append(out_f)
                        
                    elif kk[0] in All_jackpots and kk[1] not in All_jackpots: # first of pair jackpot 
                        out_f=kk[0], kk[1], k[1], m[3],  m[1], 'M', 'F', Jackpot, Non_jackpot, kkk[0], kkk[1]
                        READ_result_merged_jack_sex.append(out_f)
                        
                    elif kk[0] not in All_jackpots and kk[1] in All_jackpots: # first of pair non-jackpot
                        out_f=kk[0], kk[1], k[1], m[3],  m[1], 'M', 'F', Non_jackpot, Jackpot, kkk[0], kkk[1]
                        READ_result_merged_jack_sex.append(out_f)
                        
                    elif kk[0] not in All_jackpots and kk[1] not in All_jackpots:
                        out_f=kk[0], kk[1], k[1], m[3],  m[1], 'M', 'F', Non_jackpot, Non_jackpot, kkk[0], kkk[1]
                        READ_result_merged_jack_sex.append(out_f)
                       
                if index_male[0] == 1: #first is female and second is male 
                #print(index_male, kk)
                    if all(e in All_jackpots for e in kk): # both samples are jackpots 
                        out_f=kk[0], kk[1], k[1], m[3],  m[1], 'F', 'M', Jackpot, Jackpot, kkk[0], kkk[1]
                        READ_result_merged_jack_sex.append(out_f)
                        
                    elif kk[0] in All_jackpots and kk[1] not in All_jackpots: # first of pair jackpot 
                        out_f=kk[0], kk[1], k[1], m[3],  m[1], 'F', 'M', Jackpot, Non_jackpot, kkk[0], kkk[1]
                        READ_result_merged_jack_sex.append(out_f)
                        
                    elif kk[0] not in All_jackpots and kk[1] in All_jackpots: # first of pair non-jackpot
                        out_f=kk[0], kk[1], k[1], m[3],  m[1], 'F', 'M', Non_jackpot, Jackpot, kkk[0], kkk[1]
                        READ_result_merged_jack_sex.append(out_f)
                        
                    elif kk[0] not in All_jackpots and kk[1] not in All_jackpots:
                        out_f=kk[0], kk[1], k[1], m[3],  m[1], 'F', 'M', Non_jackpot, Non_jackpot, kkk[0], kkk[1]
                        READ_result_merged_jack_sex.append(out_f)
                        
                        
            else: #both samples are males 
                if all(e in All_jackpots for e in kk): # both samples are jackpots 
                    out_f=kk[0], kk[1], k[1], m[3],  m[1], 'M', 'M', Jackpot, Jackpot, kkk[0], kkk[1]
                    READ_result_merged_jack_sex.append(out_f)
                elif kk[0] in All_jackpots and kk[1] not in All_jackpots: # first of pair jackpot 
                    out_f=kk[0], kk[1], k[1], m[3],  m[1], 'M', 'M', Jackpot, Non_jackpot, kkk[0], kkk[1]
                    READ_result_merged_jack_sex.append(out_f)

                elif kk[0] not in All_jackpots and kk[1] in All_jackpots: # first of pair non-jackpot
                    out_f=kk[0], kk[1], k[1], m[3],  m[1], 'M', 'M', Non_jackpot, Jackpot, kkk[0], kkk[1]
                    READ_result_merged_jack_sex.append(out_f)
                elif kk[0] not in All_jackpots and kk[1] not in All_jackpots: 
                    out_f=kk[0], kk[1], k[1], m[3],  m[1], 'M', 'M', Non_jackpot, Non_jackpot, kkk[0], kkk[1]
                    READ_result_merged_jack_sex.append(out_f)

                    

        elif any(e in Males_all_samps for e in kk) == False: #both females 
            
            if all(e in All_jackpots for e in kk): # both samples are jackpots 
                out_f=kk[0], kk[1], k[1], m[3],  m[1], 'F', 'F', Jackpot, Jackpot, kkk[0], kkk[1]
                READ_result_merged_jack_sex.append(out_f)
                
            elif kk[0] in All_jackpots and kk[1] not in All_jackpots: # first of pair jackpot 
                out_f=kk[0], kk[1], k[1], m[3],  m[1], 'F', 'F', Jackpot, Non_jackpot, kkk[0], kkk[1]
                READ_result_merged_jack_sex.append(out_f)
                
            elif kk[0] not in All_jackpots and kk[1] in All_jackpots: # first of pair non-jackpot
                out_f=kk[0], kk[1], k[1], m[3],  m[1], 'F', 'F', Non_jackpot, Jackpot, kkk[0], kkk[1]
                READ_result_merged_jack_sex.append(out_f)
                
            elif kk[0] not in All_jackpots and kk[1] not in All_jackpots:
                out_f=kk[0], kk[1], k[1], m[3],  m[1], 'F', 'F', Non_jackpot, Non_jackpot, kkk[0], kkk[1]
                READ_result_merged_jack_sex.append(out_f)
        else:
            print(len(index_male),'how')
               
    else: # k[1] !='First Degree' or k[1] !='Second Degree':    
        if all(e in Males_all_samps for e in kk):  #both samples are males 
            if all(e in All_jackpots for e in kk): # both samples are jackpots 
                out_f=kk[0], kk[1], k[1], m[3],  m[1], 'M', 'M', Jackpot, Jackpot, kkk[0], kkk[1]
                
                READ_result_merged_jack_sex.append(out_f)
            elif kk[0] in All_jackpots and kk[1] not in All_jackpots: # first of pair jackpot 
                out_f=kk[0], kk[1], k[1], m[3],  m[1], 'M', 'M', Jackpot, Non_jackpot, kkk[0], kkk[1]
                READ_result_merged_jack_sex.append(out_f)
                
            elif kk[0] not in All_jackpots and kk[1] in All_jackpots: # first of pair non-jackpot
                out_f=kk[0], kk[1], k[1], m[3],  m[1], 'M', 'M', Non_jackpot, Jackpot, kkk[0], kkk[1]
                READ_result_merged_jack_sex.append(out_f)
            elif kk[0] not in All_jackpots and kk[1] not in All_jackpots:
                out_f=kk[0], kk[1], k[1], m[3],  m[1], 'M', 'M', Non_jackpot, Non_jackpot, kkk[0], kkk[1]
                READ_result_merged_jack_sex.append(out_f)
                
        if any(e in Males_all_samps for e in kk):
            index_male=[iii for iii, e in enumerate(kk) if e in Males_all_samps]
            if len(index_male) == 1:
                if index_male[0] == 0: #first is male and second is female 
                    
                    if all(e in All_jackpots for e in kk): # both samples are jackpots 
                        out_f=kk[0], kk[1], k[1], m[3],  m[1], 'M', 'F', Jackpot, Jackpot, kkk[0], kkk[1]
                        READ_result_merged_jack_sex.append(out_f)
                        
                    elif kk[0] in All_jackpots and kk[1] not in All_jackpots: # first of pair jackpot 
                        out_f=kk[0], kk[1], k[1], m[3],  m[1], 'M', 'F', Jackpot, Non_jackpot, kkk[0], kkk[1]
                        READ_result_merged_jack_sex.append(out_f)
                        
                    elif kk[0] not in All_jackpots and kk[1] in All_jackpots: # first of pair non-jackpot
                        out_f=kk[0], kk[1], k[1], m[3],  m[1], 'M', 'F', Non_jackpot, Jackpot, kkk[0], kkk[1]
                        READ_result_merged_jack_sex.append(out_f)
                        
                    elif kk[0] not in All_jackpots and kk[1] not in All_jackpots:
                        out_f=kk[0], kk[1], k[1], m[3],  m[1], 'M', 'F', Non_jackpot, Non_jackpot, kkk[0], kkk[1]
                        READ_result_merged_jack_sex.append(out_f)
                       
                if index_male[0] == 1: #first is female and second is male 
                #print(index_male, kk)
                    if all(e in All_jackpots for e in kk): # both samples are jackpots 
                        out_f=kk[0], kk[1], k[1], m[3],  m[1], 'F', 'M', Jackpot, Jackpot, kkk[0], kkk[1]
                        READ_result_merged_jack_sex.append(out_f)
                        
                    elif kk[0] in All_jackpots and kk[1] not in All_jackpots: # first of pair jackpot 
                        out_f=kk[0], kk[1], k[1], m[3],  m[1], 'F', 'M', Jackpot, Non_jackpot, kkk[0], kkk[1]
                        READ_result_merged_jack_sex.append(out_f)
                        
                    elif kk[0] not in All_jackpots and kk[1] in All_jackpots: # first of pair non-jackpot
                        out_f=kk[0], kk[1], k[1], m[3],  m[1], 'F', 'M', Non_jackpot, Jackpot, kkk[0], kkk[1]
                        READ_result_merged_jack_sex.append(out_f)
                        
                    elif kk[0] not in All_jackpots and kk[1] not in All_jackpots:
                        out_f=kk[0], kk[1], k[1], m[3],  m[1], 'F', 'M', Non_jackpot, Non_jackpot, kkk[0], kkk[1]
                        READ_result_merged_jack_sex.append(out_f)
                       

        if any(e in Males_all_samps for e in kk) == False:
            
            if all(e in All_jackpots for e in kk): # both samples are jackpots 
                out_f=kk[0], kk[1], k[1], m[3],  m[1], 'F', 'F', Jackpot, Jackpot, kkk[0], kkk[1]
                READ_result_merged_jack_sex.append(out_f)
                
            elif kk[0] in All_jackpots and kk[1] not in All_jackpots: # first of pair jackpot 
                out_f=kk[0], kk[1], k[1], m[3],  m[1], 'F', 'F', Jackpot, Non_jackpot, kkk[0], kkk[1]
                READ_result_merged_jack_sex.append(out_f)
                
            elif kk[0] not in All_jackpots and kk[1] in All_jackpots: # first of pair non-jackpot
                out_f=kk[0], kk[1], k[1], m[3],  m[1], 'F', 'F', Non_jackpot, Jackpot, kkk[0], kkk[1]
                READ_result_merged_jack_sex.append(out_f)
                
            elif kk[0] not in All_jackpots and kk[1] not in All_jackpots:
                out_f=kk[0], kk[1], k[1], m[3],  m[1], 'F', 'F', Non_jackpot, Non_jackpot, kkk[0], kkk[1]
                READ_result_merged_jack_sex.append(out_f)
        
            



file_path = f'/gpfs/scratch/akwakye/AIM_1/lcmlkin/READ/per_timepoint/READ_multi_timepoints_result_merged_jack_sex.txt'
with open(file_path, 'w') as file:
    file.write("Pair_1," +"Pair_2,"+ "Relationship,"+ "PO,"+ "Normalized_PO,"+ "Pair_1_Sex,"+ "Pair_2_Sex,"+ "Pair_1_Dosage,"+ "Pair_2_Dosage,"+ "Pair_1_Year,"+ "Pair_2_Year,"+'\n')
    for i in range(len(READ_result_merged_jack_sex)):
        line =','.join(map(str, READ_result_merged_jack_sex[i]))  # Convert each element to string and join with space
        
        file.write(line+'\n')



#the above script will produce a file that is fed into the following to plot the networks.


List_of_files=open('READv2_all_timepoints_result_with_SC2020_merged_jack_sex.txt','r') # lcMLkin_result_merged_jack_sex.txt lcMLkin_neutral_chroms_result_merged_jack_sex.txt
File_list=List_of_files.read()
READ_result=File_list.split('\n')
if READ_result[-1] =='':
    del(READ_result[-1])




G_13_14=nx.Graph()
G_all_READ=nx.Graph()
G_13_14.name='SC2013_SC2014'
READ_r_coef=[]
READ_all_connections=[]
sc_13_14_count_1=[]
sc_13_14_count_2=[]
sc_13_14_count_3=[]

for i in range(1,len(READ_result)):
    k=READ_result[i].split(',')
    key=[k[0],k[1]]
    g=int(k[-2])#k[0].split('-')[1]
    gg=int(k[-1]) #k[1].split('-')[1]
        
        
    #if k[2] == 'First Degree' or k[2] == 'Second Degree' or k[2] == 'Third Degree'  :
        
    if (int(g) == 2013 and int(gg)== 14) or (int(g) == 2013 and int(gg)== 2013) or (int(g) == 14 and int(gg)== 14) or (int(g) == 14 and int(gg)== 2013):
    #if (int(g) == 2013 and int(gg)== 14) or (int(g) == 14 and int(gg)== 2013):
    #if (int(g) == 2013 and int(gg)== 2013):
    #if (int(g) == 14 and int(gg)== 14):
        if k[2] == 'First Degree':
            degree=1
            #print(k)
            sc_13_14_count_1.append(k)
        if k[2] == 'Second Degree':
            degree=0.5
            #print(k)
            sc_13_14_count_2.append(k)
        if k[2] == 'Third Degree':
            degree=0.1
            #print(k)
            sc_13_14_count_3.append(k)
        if k[7]== 'Jackpot' and k[8] =='Jackpot':
            
            G_13_14.add_edge(key[0], key[1], weight=float(degree), category='A')
            G_13_14.nodes[k[0]]['year']=int(g)
            G_13_14.nodes[k[1]]['year']=int(gg)
            
            G_13_14.nodes[k[0]]['mito-hap']=mito_dict[key[0]]
            G_13_14.nodes[k[1]]['mito-hap']=mito_dict[key[1]]
            
            #print(k[5], k[6])
            if key[0] in NRY_haps_samps:
                G_13_14.nodes[k[0]]['y-hap']=NRY_dict[key[0]]
                
            if key[1] in NRY_haps_samps:
                G_13_14.nodes[k[1]]['y-hap']=NRY_dict[key[1]]
                
            G_13_14.nodes[k[0]]['dosage']=k[7]
            G_13_14.nodes[k[1]]['dosage']=k[8]
            G_13_14.nodes[k[0]]['category']='A'
            G_13_14.nodes[k[1]]['category']='A'
        elif (k[7]== 'Jackpot'  and k[8] =='Non_jackpot') or (k[7]== 'Non_jackpot' and k[8] =='Jackpot'): #k[5]== 'Non_jackpot' 'Non_jackpot' and k[6] =='Jackpot':
            G_13_14.add_edge(k[0],k[1], weight=float(degree), category='B')
            G_13_14.nodes[k[0]]['year']=int(g)
            G_13_14.nodes[k[1]]['year']=int(gg)
            
            G_13_14.nodes[k[0]]['mito-hap']=mito_dict[key[0]]
            G_13_14.nodes[k[1]]['mito-hap']=mito_dict[key[1]]
           
            #print(k[5], k[6])
            if key[0] in NRY_haps_samps:
                G_13_14.nodes[k[0]]['y-hap']=NRY_dict[key[0]]
                
            if key[1] in NRY_haps_samps:
                G_13_14.nodes[k[1]]['y-hap']=NRY_dict[key[1]]
           
            G_13_14.nodes[k[0]]['dosage']=k[7]
            G_13_14.nodes[k[1]]['dosage']=k[8]
            G_13_14.nodes[k[0]]['category']='B'
            G_13_14.nodes[k[1]]['category']='B'

        elif k[7]== 'Non_jackpot' and k[8] =='Non_jackpot':

            G_13_14.add_edge(k[0],k[1], weight=float(degree), category='C')
            G_13_14.nodes[k[0]]['year']=int(g)
            G_13_14.nodes[k[1]]['year']=int(gg)
            
            G_13_14.nodes[k[0]]['mito-hap']=mito_dict[key[0]]
            G_13_14.nodes[k[1]]['mito-hap']=mito_dict[key[1]]
           
            if key[0] in NRY_haps_samps:
                G_13_14.nodes[k[0]]['y-hap']=NRY_dict[key[0]]
                
            if key[1] in NRY_haps_samps:
                G_13_14.nodes[k[1]]['y-hap']=NRY_dict[key[1]]
           
            
            G_13_14.nodes[k[0]]['category']='C'
            G_13_14.nodes[k[1]]['category']='C'
            G_13_14.nodes[k[0]]['dosage']=k[7]
            G_13_14.nodes[k[1]]['dosage']=k[8]





pos = graphviz_layout(G_13_14)  
# nodes = list(pos.keys())
# positions = list(pos.values())
# random.shuffle(positions)

# # Reassign the shuffled positions to the nodes
# pos = {node: positions[i] for i, node in enumerate(nodes)}


plt.figure(figsize=(10, 10))

cat_A = [(u, v) for (u, v, d) in G_13_14.edges(data=True) if d["category"] == 'A']
cat_B = [(u, v) for (u, v, d) in G_13_14.edges(data=True) if d["category"] == 'B']
cat_C = [(u, v) for (u, v, d) in G_13_14.edges(data=True) if d["category"] == 'C']

#node_colors= {14:'#EDEA18' } #2013: '#A93C08' ,'GGC':'#AF00DB','GGGGGCG':'#C70EC8'#,2013: '#A93C08' 14:'#EDEA18'  #14:'#EDEA18' 2013: '#A93C08'
#node_colors= {2013: '#A93C08',14:'#EDEA18' } #2013: '#A93C08' ,'GGC':'#AF00DB','GGGGGCG':'#C70EC8'#,2013: '#A93C08' 14:'#EDEA18'  #14:'#EDEA18' 2013: '#A93C08'
node_colors= {2013: '#A93C08', 14:'#EDEA18' } # #2013: '#A93C08'

k=list(G_13_14.degree())
k_sort=sorted(k, key=lambda x: x[1], reverse=True)

node_types = {'Jackpot': 'o', 'Non_jackpot': '^'}
#node_types = {'Jackpot': 'o', 'Non_jackpot': '^', 'GGGGGCG':'D','GGC':'h' }

for node, attrs in G_13_14.nodes(data=True):
    year = attrs['year']
    deg=G_13_14.degree(node)
    dos=attrs['dosage']
#     if len(attrs) > 4:
#         if attrs['y-hap']=='GGC':
#             print('k', attrs)
#             nx.draw_networkx_nodes(G_13_14, pos, nodelist=[node],node_color=node_colors['GGC'] ,  node_shape=node_types['GGC'], node_size=deg*30)
            
#         elif attrs['mito-hap']=='GGGGGCG':
#             nx.draw_networkx_nodes(G_13_14, pos, nodelist=[node],node_color=node_colors['GGGGGCG'] ,  node_shape=node_types['GGGGGCG'], node_size=deg*30)
#         else:
#             nx.draw_networkx_nodes(G_13_14, pos, nodelist=[node],node_color=node_colors[year] ,  node_shape=node_types[dos], node_size=deg*30)
#             #print(attrs)
#     else:
#         nx.draw_networkx_nodes(G_13_14, pos, nodelist=[node],node_color=node_colors[year] ,  node_shape=node_types[dos], node_size=deg*30)
    
        
    nx.draw_networkx_nodes(G_13_14, pos, nodelist=[node],node_color=node_colors[year] ,  node_shape=node_types[dos], node_size=deg*30)
    
#     nx.draw_networkx_nodes(G_13_14, pos, nodelist=[node],node_color=node_colors[year] ,  node_shape=node_types[dos], node_size=deg*30)

#nx.draw_networkx_nodes(G_13_14, pos, node_size=10)

nx.draw_networkx_edges(G_13_14, pos, edgelist=cat_A, width=[d['weight']*5 for (u, v, d) in G_13_14.edges(data=True) if (u, v) in cat_A],alpha=0.8, edge_color="blue",style="solid", label='Jackpot-Jackpot')
nx.draw_networkx_edges(G_13_14, pos, edgelist=cat_B,width=[d['weight']*5 for (u, v, d) in G_13_14.edges(data=True) if (u, v) in cat_B],  alpha=0.8, edge_color="green", style="solid", label='Non_jackpot-Jackpot')
nx.draw_networkx_edges(G_13_14, pos, edgelist=cat_C,width=[d['weight']*5 for (u, v, d) in G_13_14.edges(data=True) if (u, v) in cat_C],  alpha=0.8, edge_color="red", style='solid', label='Non_jackpot-Non_jackpot')

node_handles = []
node_labels = []
for node_type, marker in node_types.items():
    node_handle = plt.scatter([], [], marker=marker, color='black', label=node_type)
    node_handles.append(node_handle)
    
for year, color in node_colors.items():
    if year ==14:
        year= 2014
        
    node_handle = plt.scatter([], [], color=color, label=str(year))
    node_handles.append(node_handle)
    node_labels.append(str(year))

edge_handles, edge_labels = plt.gca().get_legend_handles_labels()

#nx.draw_networkx_labels(G_13_14, pos, labels={node: node for node in G_13_14.nodes()}, font_size=4.5, font_color='black', font_weight='bold')


#nx.draw_networkx_labels(G_13_14, pos, labels={k_sort[0][0]: k_sort[0][0] }, font_size=10, font_color='black', font_weight='bold')
#nx.draw_networkx_labels(G_13_14, pos, labels={'SC-2013-PB-X140-G5': 'SC-2013-PB-X140-G5' }, font_size=10, font_color='black', font_weight='bold')

# nx.draw_networkx_labels(G_13_14, pos, labels={k_sort[0][0]: k_sort[0][0] }, font_size=10, font_color='black', font_weight='bold')



# for node, (x, y) in pos.items():
#     plt.text(x, y, str(node), fontsize=10, ha='center', va='center', rotation=270, rotation_mode='anchor', fontweight='bold')


plt.legend(handles=node_handles, labels=node_labels,   fontsize=10)
plt.legend(handles=edge_handles, labels=edge_labels,  fontsize=10)
#plt.title('Scout 2014')
plt.title('SC2013 and SC2014',fontsize=16)
#plt.title('SC2014',fontsize=16)

plt.savefig('sc_14_READ_v2.pdf', format='pdf', dpi=600, bbox_inches='tight')
#plt.savefig('sc_13_14_READ_v2.pdf', format='pdf', dpi=600, bbox_inches='tight')

#plt.savefig('sc_13_READ_v2.pdf', format='pdf', dpi=600, bbox_inches='tight')

plt.show()
plt.close()







