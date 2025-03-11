
import time
import vcf 
import math 
import numpy as np
import pandas as pd
from collections import Counter 
import matplotlib.pyplot as plt
import numpy.ma as ma
import matplotlib.pyplot as plt
import os
import pysam
import copy


rho_cut_off=0.8
dist_cut_off=1000
rho_cum_cut_off=0.80
rho_strict=0.99
min_freq=0.1
max_freq=1-min_freq

filt_missingness=2.62 #this means a site that is present in at least 1/{value} of the samples is included. 

t0 = time.time()
All_males_vcf_reader = vcf.Reader(open('All_males_RS_SC_13_14_15_17_gatk_bed_list_interval_ploidy_1.removeindel_pad200.vcf.gz', 'r'))

nucleotides=["A","C", "G","T", None]
All_samps=[]
positions=[]
SNP_info=[]
for record in All_males_vcf_reader:
    SNP_info.append([record.INFO["AF"], record.INFO["AC"], record.INFO["AN"]])
    #print(record.INFO["AN"],record.INFO["AC"], len(record.samples)/filt_missingness)
    if len(record.INFO["AF"]) == 1 and record.INFO["AN"] > (len(record.samples)/filt_missingness): #filtering for common variants 
        Sample_genotypes=record.samples
        sample_genos=[]
        samp_names=[]
        for call in Sample_genotypes:
            base = call.gt_bases
            #print(base)
            pos=call.site.POS
            samples=call.sample
            sample_genos.append(base)
            samp_names.append(samples)
        All_samps.append(sample_genos) # all sites and genotypes of all samples
        positions.append(pos)
    

t1 = time.time()
print('finished loading sites', len(All_samps), "time elapsed: {} s".format(t1-t0))


t0 = time.time()
All_Samp_pos=list(zip(positions, All_samps)) #all sites included with genotypes for all individuals 
All_Samp_pos_ar=np.array(All_Samp_pos, dtype=object)
All_Samp_pos_ar_transposed=np.transpose(All_Samp_pos_ar)

print(len(All_Samp_pos))

site_geno_len=[]
site_genotypes=[]
positions=[]
site_genotypes_to_use=[]
for i in range(len(All_Samp_pos)):
    site_genos=[]
    site_genos_to_use=[]
    for j in range(len(All_Samp_pos[i][1])):
        if All_Samp_pos[i][1][j]!=None:
            pos=All_Samp_pos[i][0]
            geno=All_Samp_pos[i][1][j]
            site_genos.append([j,geno])
            site_genos_to_use.append(geno)
    positions.append(pos)
    site_genotypes.append(site_genos) # this is a list of individuals with called genotypes and the genotypes for all sites
    site_genotypes_to_use.append(site_genos_to_use)
    site_geno_len.append(len(site_genos))
    
informative_snps=[] #informative snps are those sites that are biallelic, and large numbers of individuals with called genotypes
informative_snps_raw=[]
informative_snps_no_impute=[]
informative_sites_genos=[]
informative_positions=[]
Base_count_freq=[]
emp_genos=[]
emp_genotypes_only=[]
for i in range(len(site_genotypes)): #filtering for maf > 0.1
    Base_count=np.unique(site_genotypes_to_use[i], return_counts=True)

    if len(Base_count[0])==2 and len(Base_count[0][1])== 1:
        
        if Base_count[1][1] > Base_count[1][0]:
            maf_geno=Base_count[0][0], Base_count[0][1]
            maf=Base_count[1][0]/sum(Base_count[1])
            emp_genos.append([str(Base_count[0][0])+"/"+str(Base_count[0][1]), maf])
            
        elif Base_count[1][1] < Base_count[1][0]:
            maf_geno=Base_count[0][1], Base_count[0][0]
            #print('k', maf_geno)
            maf=Base_count[1][1]/sum(Base_count[1])
            
        if maf > 0.1:
            emp_genos.append([str(Base_count[0][1])+"/"+str(Base_count[0][0]), maf])
            emp_genotypes_only.append(maf_geno) # to be used for dosage calculation
            
            informative_snps.append([positions[i],site_genotypes[i]]) #contains only individuals with called genotypes
            informative_sites_genos.append(site_genotypes_to_use[i]) #contains only individual alleles
            
            informative_snps_raw.append(All_Samp_pos[i][1])  # has all individuals with called and uncalled genotypes
            informative_snps_no_impute.append(All_Samp_pos[i][1])
            informative_positions.append(positions[i])

informative_snps_no_imp_ar=np.array(informative_snps_no_impute)

informative_snps_no_impute_deep_copy= informative_snps_no_impute

all_pos=[]
for i in range(len(informative_snps)): 
    snps=informative_snps[i][1]
    for j in range(len(snps)):
        pos=snps[j][0]
        all_pos.append(pos)
    

nucleotide_per_site=[]
for i in range(len(informative_snps)):
    info=informative_snps[i][1]
    nucleos=[]
    for j in range(len(info)):
        nucleos.append(info[j][1])
    nucleotide_per_site.append(nucleos)

  
    
#calculate the rate of missingness before imputation 
sites_to_use=[]
missingness_pre_imputation_dict={}
for i in range(len(informative_snps)):
    for j in range(len(All_samps)):
        if informative_snps[i][0]==All_Samp_pos[j][0]:
            m=np.sum(np.array(All_Samp_pos[j][1]) ==None)
            m_rate=m/len(All_Samp_pos[j][1])
                        
            missingness_pre_imputation_dict[All_Samp_pos[j][0]]=m_rate, m 
    
            sites_to_use.append(All_samps[j])

            
            
            
            
            
# The following will transform the genotypes into doses. 0 for first allele and 1 for second allele 
sites_to_use_dosage=[]
for i in range(len(informative_snps_raw)):
    sites=[]
    for j in range(len(informative_snps_raw[i])):
        x=informative_snps_raw[i][j]
        #print(x)
        if x==emp_genotypes_only[i][0]: 
            #print(x)
            sites.append(0)
        elif x==emp_genotypes_only[i][1]:
             sites.append(1)
        else: 
             sites.append(np.nan)
    sites_to_use_dosage.append(sites)

sites_to_use_ar=np.array(sites_to_use_dosage)

masked_sites_to_use_ar=np.ma.masked_invalid(sites_to_use_ar)


t1 = time.time()
print('biallelic sites',len(sites_to_use), "time elapsed: {} s".format(t1-t0))


t0 = time.time()
def r_squared(A,B):
    
    mu_A=np.mean(A.astype(int))
    fil_A=A.astype(int)
    dif_A=fil_A-mu_A
    diff_A=dif_A*dif_A
    var_A=np.sum(diff_A[~np.isnan(diff_A)])/(len(diff_A[~np.isnan(diff_A)])-1)
    var_A

    mu_B=np.mean(B.astype(int))
    fil_B=B.astype(int)
    dif_B=fil_B-mu_B
    diff_B=dif_B*dif_B
    var_B=np.sum(diff_B[~np.isnan(diff_B)])/(len(diff_B[~np.isnan(diff_B)])-1)
    var_B
    diff=dif_A*dif_B
    cov_A_B=np.sum(diff[~np.isnan(diff)])/(len(diff[~np.isnan(diff)])-1)
    cov_A_B

    r=cov_A_B/(np.sqrt(var_A*var_B))
    r_squared= r ** 2
    return r



rho_mat=[]

for i in range(len(masked_sites_to_use_ar)):
    loci_rho=[]
    for j in range(len(masked_sites_to_use_ar)):
        rho=r_squared(masked_sites_to_use_ar[i],masked_sites_to_use_ar[j])
        loci_rho.append(rho)
    rho_mat.append(loci_rho)
                
rho_mat_ar=np.array(rho_mat)   


t1 = time.time()
print('finished estimating correlation matrix', "time elapsed: {} s".format(t1-t0))


t0 = time.time()
loci_filtered_0_81=[]
loci_all_r=[]
loci_r=[]
locus_for_impute=[]
for i in range(len(rho_mat_ar)):
    k=rho_mat_ar[i]
    kk=[]
    g=[]
    gg=[]
    f=[]
    info_sites=[]
    for j in range(len(k)):
        if k[j] > rho_cut_off or k[j] < -rho_cut_off:
            #print(i, k[j], j)
            g.append([informative_positions[i], k[j], informative_positions[j]])
            f.append([i, informative_positions[i], k[j], informative_positions[j], j])
            gg.append([i, k[j], j])
            #print(g)
            kk.append(k[j])
            info_sites.append(informative_positions[j])
    loci_r.append(g)
    loci_all_r.append(f)
    loci_filtered_0_81.append(kk)
    locus_for_impute.append(gg)



loci_to_use_group=[]
actual_loci=[]
loci_all_r_grp=[]
for i in range(1, len(loci_filtered_0_81)):
    if len(loci_filtered_0_81[i]) > 1:
        actual_loci.append(loci_r[i])
        loci_all_r_grp.append(loci_all_r[i])
        loci_to_use_group.append(loci_filtered_0_81[i])
       
       
filt_loci=[]
for i in range(len(locus_for_impute)):
    if len(locus_for_impute[i]) > 1: 
        #print(locus_for_impute[i])
        filt_loci.append(locus_for_impute[i])

        
                
def group_by_length(input_list):
    groups = []
    current_group = [input_list[0]]
    for i in range(1, len(input_list)):
        if len(input_list[i]) == len(input_list[i - 1]):
            current_group.append(input_list[i])
        else:
            groups.append(current_group)
            current_group = [input_list[i]]

    groups.append(current_group)

    return groups

grouped_loci=group_by_length(filt_loci)
grouped_actual_loci=group_by_length(actual_loci)

grouped_loci_filt=[]
grouped_actual_loci_filt=[]
for i in range(len(grouped_loci)):
    
    grouped_actual_loci_filt.append(grouped_actual_loci[i][0])
    grouped_loci_filt.append(grouped_loci[i][0])
    
informative_snps_no_imp_ar_t=np.transpose(informative_snps_no_imp_ar)



t1 = time.time()
print('prepare for imputation',"time elapsed: {} s".format(t1-t0))
#imputation 

t0 = time.time()
print('start imputation')



grouped_loci_new=[]
for m in range(len(grouped_loci_filt)):
    grp=grouped_loci_filt[m]
    #print(len(grp))
    new_group=[]
    for n in range(1,len(grp)):
        new_group.append(grp[n])
        
        #if grp[n][0] !=grp[n][2]:
            #print(grp[n]) 
         #   new_group.append(grp[n])
    grouped_loci_new.append(new_group)
    #group=grouped_loci_filt[m] # take one group of sites 

#extract genotype before imputation 
lll=[]
for m in range(len(grouped_loci_new)):
    group=grouped_loci_new[m] # take one group of sites 
    group_sort=sorted(group, key=lambda x: x[1], reverse=True)  #sort them according to the cor
    for l in range(len(group_sort)):
        loci=group_sort[l]
        #print(len(group_sort))
        lll.append(loci)
        #ind_genotypes=[]
        
        
        
loci_to_extract=[]
for i in range(len(lll)):
    loci_to_extract.append(lll[i][2])
    
loci_to_ext=np.unique(loci_to_extract)

pop_geno_pre_imp=[]
for i in range(len(informative_snps_no_imp_ar_t)):
    samp_genos_pre_imp=[]
    samp_genos=informative_snps_no_imp_ar_t[i]
    for j in range(len(loci_to_ext)):
        samp_genos_pre_imp.append(samp_genos[loci_to_ext[j]])
    pop_geno_pre_imp.append(samp_genos_pre_imp)
    
### start imputation #### 
sites_to_use_dosage_t=np.transpose(sites_to_use_dosage)
informative_snps_no_imp_ar_t_dpcp= np.copy(informative_snps_no_imp_ar_t)

lll=[]
cons_all=[]
nucleos=["A","C", "G","T"]
for a in range(len(informative_snps_no_imp_ar_t)):
    gen_dict={}
    ind=informative_snps_no_imp_ar_t[a]
    ind_to_be_imputed=informative_snps_no_imp_ar_t[a]
    
    for m in range(len(grouped_loci_new)):
        group=grouped_loci_new[m]
        hggg=[]
        All_samp_imputed_geno_2=[]
        group_sort=sorted(group, key=lambda x: abs(x[1]), reverse=True)  #sort them according to the cor
        cons_nuc=[]
        loci_lll=[]
        rho_lll=[]
        dosage_lll=[]
        emp_geno_lll=[]
        pos_lll=[]
        for l in range(len(group)):
            loci=group[l]
            loci_lll.append(loci[2])
            dosage_lll.append(sites_to_use_dosage_t[a][loci[2]])
            rho_lll.append(rho_mat_ar[loci[0]][loci[2]])
            cons_nuc.append(ind[[loci[2]]][0])
            emp_geno_lll.append(emp_genotypes_only[loci[2]])
            pos_to_imp=informative_positions[loci[0]]
            pos_lll.append(informative_positions[loci[2]])
            #print(cons_nuc)
            if ind[[loci[0]]]== None:
                if any(c is not None for c in cons_nuc): #if all the sites in the group are none, they are useless 
                    dist=[abs(e - pos_to_imp) for e in pos_lll] # calculates the physical distance between the site to be predicted and the predicting sites 
                    if any(e > dist_cut_off for e in dist): # filter for physical distance 
                        indices_dist=[iii for iii, e in enumerate(dist) if e > dist_cut_off] # make this a parameter 
                        
                        cons_nuc_dist_filt = [cons_nuc[index] for index in indices_dist]
                        rho_dist_filt=[rho_lll[index] for index in indices_dist]
                        emp_dist_filt= [emp_geno_lll[index] for index in indices_dist]
                        dosage_dist_filt=[dosage_lll[index] for index in indices_dist]
                        #print(cons_nuc_dist_filt, indices_dist)
                        
                        if len(cons_nuc_dist_filt) > 1: # then filter for whether there are multiple predictors.
                            indices=[index for index, e in enumerate(cons_nuc_dist_filt) if e != None] # if there are multiple, which ones are called and which ones aren't 
                            
                            cons_nuc_indices=[cons_nuc_dist_filt[index] for index in indices] # extract the called nucleotides
                            rho_indices=[rho_dist_filt[index] for index in indices] #extract the correlation coefficient for those sites 
                            emp_geno_indices= [emp_dist_filt[index] for index in indices] # extract the minor and major alleles
                            dosage_indices=[dosage_dist_filt[index] for index in indices] # extract the dosgae of the minor and major alleles (remember, 0 is minor, 1 is major)
                            dosage_indices_rep=[-1.0 if e == 0.0 else e for e in dosage_indices] # replace zeros with -1 
                            rho_indices_abs=[abs(rho_dist_filt[index]) for index in indices]
                            if len(rho_indices)> 0:
                                rho_dos_mult= [a * b for a, b in zip(rho_indices_abs, dosage_indices_rep)] # calculates the weights of the correlations and dosgae 
                                rho_dos_mult_average=sum(rho_dos_mult)/len(rho_dos_mult)
                                
                            if all(e == cons_nuc_indices[0] for e in cons_nuc_indices): #check if all elements in a list are the same(whether all sites predict the same nucleotide) 
                                
                                if rho_dos_mult_average > rho_cum_cut_off and len(rho_indices)> 0: # if all sites are all positively correlated to the site on average, its an easy substitution

                                    ind[[loci[0]]]=cons_nuc_indices[0]
                                    print('imputed')
                                    #print(ind[[loci[0]]], cons_nuc_indices[0])
                                elif rho_dos_mult_average < -rho_cum_cut_off and len(rho_indices)> 0: # they are negatively correlated on average

                                    if all(e == emp_geno_indices[0] for e in emp_geno_indices): #checks if the same alleles are at all predicting sites
                                        
                                        if emp_geno_indices[0][0] == np.unique(cons_nuc_indices)[0]:
                                            al_to_use=emp_geno_indices[0][1]
                                            ind[[loci[0]]]=al_to_use 
                                        else:
                                            #print(emp_geno_indices[0][0])
                                            'if this prints, algorithm isnt working'
                                    else:
                                        'kkkk nothung here '
                                        #print(emp_geno_indices, '||', cons_nuc_indices)
                                        
                                
                            else: # if the multiple predicting sites disagree (think about including this later)
                                if rho_dos_mult_average > rho_cum_cut_off and len(rho_indices)> 0:
                                    if rho_indices[0] > rho_strict: # use the site with the highest correlation if its almost perfect correlation . 
                                        al_to_use=cons_nuc_indices[0]
                                        ind[[loci[0]]]=al_to_use
                                        'kkk'
                                elif rho_dos_mult_average < -rho_cum_cut_off and len(rho_indices)> 0:
                                    if rho_indices[0] < -rho_strict:
                                        if emp_geno_indices[0][0] == np.unique(cons_nuc_indices)[0]:
                                            al_to_use=emp_geno_indices[0][1]
                                            ind[[loci[0]]]=al_to_use
                                        #al_to_use=cons_nuc_indices[0]
                                            'kkk'
                               
                        else:
                            #print(cons_nuc_dist_filt)
                            if cons_nuc_dist_filt[0] != None and rho_dist_filt[0] > rho_cut_off and len(rho_indices)> 0: #imputation if there's only one prediciting site and is positively correlated 
                                #print(cons_nuc_indices[0], emp_geno_indices[0],  rho_indices[0])
                                
                                ind[[loci[0]]]=cons_nuc_dist_filt[0]
                                #print(cons_nuc_dist_filt[0])
                                'kkk'
                                
                            elif cons_nuc_dist_filt[0] != None and rho_dist_filt[0] < -rho_cut_off and len(rho_indices)> 0:
                                'kkk' 
                                if rho_dist_filt[0] < -rho_strict:
                                    if emp_dist_filt[0][0] == cons_nuc_dist_filt[0]:
                                        al_to_use=emp_geno_indices[0][1]
                                        ind[[loci[0]]]=al_to_use
                                        #print(cons_nuc_dist_filt[0], rho_dist_filt, emp_dist_filt)
                                        'kkk' 
                                        
                                        


informative_snps_post_impute=np.transpose(informative_snps_no_imp_ar_t)
informative_snps_pre_impute=np.transpose(informative_snps_no_imp_ar_t_dpcp)
sites_with_min=[]
for i in range(len(informative_snps_post_impute)):
    
    if np.sum(informative_snps_post_impute[i]==None)< 50:
        print(np.sum(informative_snps_post_impute[i]==None))
        sites_with_min.append(informative_snps_post_impute[i])
sites_with_min_t=np.transpose(sites_with_min)


for_phylo_tree=[]
all_males_=[]
selected_samples=[]
for i in range(len(sites_with_min_t)):
    if all(sites_with_min_t[i]!=None):
        #print(samp_names[i],sites_with_min_t[i])
        k=samp_names[i],sites_with_min_t[i]
        for_phylo_tree.append(k)
        all_males_.append(sites_with_min_t[i])
        selected_samples.append(samp_names[i])

all_males_samp_names =selected_samples

all_males_samples_to_use=[]
Scouts_t=(np.array(all_males_))
file_out=f"all_males_to_use_{filt_missingness}_{len(sites_with_min)}_sites"
with open(file_out, 'w') as file:
    for i in range(len(Scouts_t)):
        S=Scouts_t[i].tolist()
        T=str(S).replace("'", '').replace(',', '')
        Z=T.replace("None", "n")
        
        u=Z.replace(" ", "")
        #print(u)
        if u.count('n') < 1: #for each site, it is allowed to be missing in x individuals
            #print(u)
            k=['' if x == '['  else x for x in u]
            kk=k[:-1]
            kkk=[x for x in kk if x !='']
            all_males_samples_to_use.append([all_males_samp_names[i],u])
            T=str(kkk).replace("'", '').replace(',', '')
            v=T.replace(" ", " ")
            v=all_males_samp_names[i]+'_'+v
            file.write(v+'\n')

            
            
            
            

all_males_fasta=all_males_samples_to_use
all_males_seqs=[]
for i in range(len(all_males_fasta)):
    a=all_males_fasta[i][1]
    k=['' if x == '[' or x == ']'  else x for x in a]
    del (k[0])
    del (k[-1])
    f=''.join(k)
    #print([Scout_13_14[i][0],f])
    all_males_seqs.append([all_males_fasta[i][0],f])



def write_multifasta(file_path, sequence_list):
    with open(file_path, 'w') as fasta_file:
        for entry in sequence_list:
            header = f">{entry[0]}\n"
            sequence = entry[1] + "\n"
            #print(len(sequence))
            fasta_file.write(header)
            fasta_file.write(sequence)


output_file_path = f'fasta_post_imputation_{len(all_males_seqs)}_{filt_missingness}_{len(sites_with_min)}_sites.fasta'

write_multifasta(output_file_path, all_males_seqs)


file_path = f'informative_sites_{filt_missingness}_{len(sites_with_min)}.txt'
with open(file_path, 'w') as file:
    for i in range(len(sites_with_min_t)):
        
        line =samp_names[i]+' '+' '.join(map(str, sites_with_min_t[i]))  # Convert each element to string and join with space
        file.write(line + '\n')
