# this notebook adds the sex, jackpot status to each sample 


import numpy as np 

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


def normalize_key(key):
    return tuple(sorted(key))
    
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








