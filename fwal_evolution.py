import numpy as np 
import matplotlib.pyplot as plt

snp_file = open('/gpfs/scratch/akwakye/jupyter_works/fw_haps_poolSNPs.CH_SC_LB.2500_50.50000.phy_specific_BF.snpnb1.res', 'r') # this file can be downloaded from https://figshare.com/s/d7e8318c2ca1cd293a13 or DOI 10.6084/m9.figshare.28566707

snp_data=snp_file.read()
snp_data=snp_data.split('\n')
snp_file.close()
if snp_data[-1]=='':
    del(snp_data[-1])
    
tempo_peaks=[]
for i in range(1,len(snp_data)):
    
    k=snp_data[i].split('\t')
    tempo_peaks.append([k[0], k[1]])
    
tempo_peaks_ar=np.array(tempo_peaks)   

# The following chunk loads the data. The data consists of an n by m array with n being the number of individuals in the timepoint and m=341, the total number of genome wide. See supplementary information section 6: DETERMINING STATES OF LOCI WITH FRESHWATER-ADAPTIVE ALLELES FROM LOW-COVERAGE SEQUENCE DATA for details on how these loci were generated. 

# The numpy arrays can be found here https://figshare.com/s/723da552af1acb3a320c or DOI 10.6084/m9.figshare.28566623

RS_19_ma=np.load('RS_19_ma_count.npy')
RS_19_samps=np.load('RS_19_Samp_names.npy')
RS_19_dict = dict(zip(RS_19_samps, RS_19_ma))

SC_13_ma=np.load('SC_13_ma_count.npy')
SC_13_samps=np.load('SC_13_Samp_names.npy')
SC_13_dict = dict(zip(SC_13_samps, SC_13_ma))

SC_14_ma=np.load('SC_14_ma_count.npy')
SC_14_samps=np.load('SC_14_samp_names.npy')
SC_14_dict = dict(zip(SC_14_samps, SC_14_ma))

SC_15_ma=np.load('SC_15_ma_count.npy')
SC_15_samps=np.load('SC_15_Samp_names.npy')
SC_15_dict = dict(zip(SC_15_samps, SC_15_ma))


SC_17_ma=np.load('SC_17_ma_count.npy')
SC_17_samps=np.load('SC_17_Samp_names.npy')
SC_17_dict = dict(zip(SC_17_samps, SC_17_ma))

SC_20_ma=np.load('SC_20_ma_count.npy')
SC_20_samps=np.load('SC_20_Samp_names.npy')
SC_20_dict = dict(zip(SC_20_samps, SC_20_ma))


# defining functions to be used for this analyses 

def count_sequences(arr, loci):
    ''' 
    This function counts the number of contiguous loci in each individual. 
    ''''
    sequences = []
    current_sequence = []
    ecopeak_contiguous=[]
    current_ecopeak=[]
    for i in range(len(arr)):
        if arr[i] == 1 or arr[i] == 0:
            chromo=loci[i].split('\t')[0]
            k=loci[i].split('\t')[3].split(',')
            kk=chromo,int(k[0]), int(k[-1])
            current_ecopeak.append(kk)
            current_sequence.append(arr[i])
        else:
            if current_sequence: 
                sequences.append(current_sequence)
                ecopeak_contiguous.append(current_ecopeak) 
                #print(sequences, ecopeak_contiguous)
#                 k=loci[i].split('\t')[3].split(',')
#                 kk=k[0], k[-1]
#                 ecopeak_contiguous.append(kk)
                current_sequence = []
                current_ecopeak=[]
    if current_sequence:
        sequences.append(current_sequence)
        ecopeak_contiguous.append(current_ecopeak)
    return sequences, ecopeak_contiguous

 
def cM_Mb(ma_alleles, snp_data, samps_name, n):

     ''' 
    This function was used to generate supplementary table 2. The output from this function was also used to plot figure 2 and supplementary figure 3. The files required for this function includes the numpy arrays containing the genotypes of each individual at all 341 loci. These files have been deposited in figshare and can be accessed with the following link: https://figshare.com/s/723da552af1acb3a320c  
    
    The map files that contain the recombination maps generated from the ancestral Rabbit Slough population can be found through this link:https://figshare.com/s/a8798163a29fc3b7a10d or with DOI 10.6084/m9.figshare.28566575 
    
    ''''

    file=open('jackpots.txt','r') # lcMLkin_result_merged_jack_sex.txt lcMLkin_neutral_chroms_result_merged_jack_sex.txt
    File_list=file.read()
    jackpots_list=File_list.split('\n')
    if jackpots_list[-1] =='':
        del(jackpots_list[-1])

    mask=ma_alleles!=-9
    ma_alleles_count_new=ma_alleles[mask].reshape(ma_alleles.shape[0], -1)

    indices=np.where(ma_alleles[0] != -9)
    tempo_peaks_filtered=np.array(snp_data[1:])[indices]

    all_sequences = []
    jackpots_l=[]
    non_jackpots_l=[]
    for idx, arr in enumerate(ma_alleles_count_new):
        k, m = count_sequences(arr, tempo_peaks_filtered)
        #print(k)
        if len(k) !=0:
            for j in range(len(k)):
                if len(k[j]) > n:
                    #print(len(k), k[j], m[j])
                    all_sequences.append([idx,k[j], m[j]])
                    if samps_name[idx] in jackpots_list:
                        jackpots_l.append([idx,k[j], m[j]])
                    elif samps_name[idx] not in jackpots_list:
                        non_jackpots_l.append([idx,k[j], m[j]])

    year=samps_name[0].split('-')[1]
    filepath='SC'+year+'_contiguous_FWALs_cM_Mb.txt'
    with open(filepath, 'w') as file:
        file.write('Sample\tGenotype_of_FWALs\tChromosome\tFWALs\tSize(In cM)\tSize(In Mb)\tNumber_of_FWALs\n')

        Genetic_physical_distance=[]
        jackpots={}
        non_jackpots={}
        jackpots_cM={}
        non_jackpots_cM={}
        Genetic_physical_distance_jacks=[]
        Genetic_physical_distance_non_jacks=[]
        count_hets=[]
        for i in all_sequences:
            #if len(i[2]) ==10:
            #print(i[1])
            count_hets.append(i[1])
            kk=i[2][-1][2], i[2][0][1]
            kkk=np.sort(kk)

            chromo=i[2][0][0].split(' ')[0]

            #kkkk=i[2][-1][2]-i[2][0][2]
            if len(i[1]) ==1:
                
                kkkk=i[2][0][2]-i[2][0][1]
            elif len(i[1]) > 1:
                kkkk=i[2][-1][2]-i[2][0][2]
            
            Mb=abs(round(kkkk/1000000,5))

            if samps_name[i[0]] in jackpots_list:
                if SC_13_samps[i[0]] not in jackpots.keys():
                    jackpots[samps_name[i[0]]]=[ Mb]
                elif SC_13_samps[i[0]] in jackpots.keys():
                    jackpots[samps_name[i[0]]].append(Mb)
            elif samps_name[i[0]] not in jackpots_list:
                if SC_13_samps[i[0]] not in non_jackpots.keys():
                    non_jackpots[samps_name[i[0]]]=[ Mb]
                elif SC_13_samps[i[0]] in non_jackpots.keys():
                    non_jackpots[samps_name[i[0]]].append(Mb)
                
                
            cM_file = open('/gpfs/scratch/akwakye/AIM_1/imputation_stickleback_TN5/CentiMorgan_map/filtered_map_files/Chrom_'+chromo+'_filtered.gmap', 'r') #remember to change this to the path that contains the map files.
            cM_data=cM_file.read()
            cM_data=cM_data.split('\n')
            if cM_data[-1]=='':
                del(cM_data[-1])
            loci_cm=[]
            for j in cM_data:
                f=j.split('\t')
                if int(f[1]) > kkk[0] and int(f[1]) < kkk[1]:
                    #print(f)
                    loci_cm.append(f)
            if len(loci_cm)> 0:
                Genetic_distance=abs(round(float(loci_cm[-1][2])-float(loci_cm[0][2]),5))
            else:
                'k'
            #print(Genetic_distance,Mb,loci_cm[0], loci_cm[-1])
                if samps_name[i[0]] in jackpots_list:
                    if SC_13_samps[i[0]] not in jackpots_cM.keys():
                        jackpots_cM[samps_name[i[0]]]=[ Mb]
                    elif SC_13_samps[i[0]] in jackpots_cM.keys():
                        jackpots_cM[samps_name[i[0]]].append(Mb)
                elif samps_name[i[0]] not in jackpots_list:
                    if SC_13_samps[i[0]] not in non_jackpots_cM.keys():
                        non_jackpots_cM[samps_name[i[0]]]=[ Mb]
                    elif SC_13_samps[i[0]] in non_jackpots_cM.keys():
                        non_jackpots_cM[samps_name[i[0]]].append(Mb)
               
            Genetic_physical_distance.append([Genetic_distance, Mb, len(i[1])])
            #print(str(samps_name[i[0]]), jackpots)
            if samps_name[i[0]] in jackpots_list:
                #print(str(samps_name[i[0]]))
                Genetic_physical_distance_jacks.append([Genetic_distance, Mb, len(i[1])])
            elif samps_name[i[0]] not in jackpots_list:
                #print('kk',samps_name[i[0]])
                Genetic_physical_distance_non_jacks.append([Genetic_distance, Mb, len(i[1])])
                
            a=','.join(map(str, i[1]))
            b=','.join(map(str, i[2]))
            #print(type(Mb))
            #print(samps_name[i[0]], a,chromo,  b, Genetic_distance, Mb, len(i[1]))
            file.write(str(samps_name[i[0]])+'\t'+str(a)+'\t'+str(chromo)+'\t'+str(b)+'\t'+str(Genetic_distance)+'\t'+ str(Mb)+'\t'+str(len(i[1]))+'\n')
            cM_mB_ar=np.array(Genetic_physical_distance)
            cM_mB_ar_jacks=np.array(Genetic_physical_distance_jacks)
            cM_mB_ar_non_jacks=np.array(Genetic_physical_distance_non_jacks)
    return cM_mB_ar,cM_mB_ar_jacks,cM_mB_ar_non_jacks, count_hets #jackpots, non_jackpots, jackpots_cM, non_jackpots_cM


toplot_rs19=cM_Mb(RS_19_ma, snp_data, RS_19_samps,0)
toplot_13=cM_Mb(SC_13_ma, snp_data, SC_13_samps,0) 
toplot_14=cM_Mb(SC_14_ma, snp_data, SC_14_samps,0)
toplot_15=cM_Mb(SC_15_ma, snp_data, SC_15_samps, 0)
toplot_17=cM_Mb(SC_17_ma, snp_data, SC_17_samps, 0)
toplot_20=cM_Mb(SC_20_ma, snp_data, SC_20_samps, 0)



# all four columns version 1 with fw content first 

def plot_graphs(datasets, labels, fw_content_all_timepoints):
    col_dict = {
        'RS2019': '#9F2305',
        'SC2013': '#A93C08',
        'SC2014': '#b3550a',
        'SC2015': '#edea18',
        'SC2017': '#819d4e',
        'SC2020': "#155084"
    }
    
    fig, axs = plt.subplots(6, 4, figsize=(20, 20))  # Adjust the figsize as needed
    axs = axs.flatten()  # Flatten the grid to easily iterate through
    
    bin_range = (-0.5, 16)
    bins = np.linspace(bin_range[0], bin_range[1], 30)  # Adjust the number of bins as needed
    
    for i, (toplot, label) in enumerate(zip(datasets, labels)):
        # Plot the histogram for fw_content_all_timepoints in the first column
        ax = axs[4 * i]  # First column index
        ax.hist(fw_content_all_timepoints[i], bins=30, color=col_dict[label])
        
        # Calculate the mean and max for the current list
        mean_fw_content = np.mean(fw_content_all_timepoints[i])
        max_fw_content = np.max(fw_content_all_timepoints[i])
        
        ax.axvline(mean_fw_content, color='red', linestyle='--', linewidth=3, label=f'Mean: {mean_fw_content:.2f}')
        
        ax.set_xlim(0, 300)
        
        if i == len(datasets) - 1:
            ax.set_xlabel('FW Content', fontsize=14)
        else:
            ax.set_xlabel('')
        
        # Set y-axis label for the first column
        ax.set_ylabel(f'{label}', fontsize=14)
        
        ax.legend(loc='upper right')
    
        for j in range(len(toplot[0][0])):
            if j == 0:
                k = 'cM'
            elif j == 1:
                k = 'Mb'
            elif j == 2:
                k = 'Number_of_FWALs'
            
            means_fwals = np.mean(toplot[0][:, j])
            
            ax = axs[4 * i + j + 1]  # Adjust index to account for the shifted column
            ax.hist(toplot[0][:, j], bins=bins, color=col_dict[label])
            ax.axvline(means_fwals, color='red', linestyle='--', linewidth=3, label=f'Mean: {means_fwals:.2f}')
            ax.set_xlim(bin_range)
            
            if i == len(datasets) - 1:
                ax.set_xlabel(f'Contiguous FWAL({k})', fontsize=14)
            else:
                ax.set_xlabel('')
            
            if (4 * i + j + 1) % 4 != 0:  # Columns 2, 3, and 4
                ax.set_ylabel('')
            else:
                ax.set_ylabel(f'{label}', fontsize=14)
            
            ax.legend(loc='upper right')
    
    plt.tight_layout()

    plot_name = f'contiguous_FWALs_hist_grid_15_with_fw_content.pdf'
    plt.savefig(plot_name, dpi=600)
    plt.show()
    plt.close()


datasets = [toplot_rs19, toplot_13, toplot_14, toplot_15, toplot_17, toplot_20]
labels = ['RS2019', 'SC2013', 'SC2014', 'SC2015', 'SC2017', 'SC2020']


plot_graphs(datasets, labels, fw_content_all_timepoints)





#### the output from the function cM_Mb was also used to calculate the changes in homozygosity over time (supplementary figure 4). To run this for each time point, simply change the timepoint argument and hets line to match the respective time point. 


timepoint='rs19' # change this line 
from itertools import chain
import matplotlib.pyplot as plt

hets=list(chain(*toplot_rs19[-1])) # change this line
hets_counts=np.unique(hets, return_counts=True)
normalized_counts = hets_counts[1] / np.sum(hets_counts[1])

plt.bar(hets_counts[0], hets_counts[1]/np.sum(hets_counts[1]), color='#155084')
print(hets_counts[1]/np.sum(hets_counts[1]))
plt.xticks([0,1],['Homozygous','Heterozygous'],  fontfamily='Sans serif', fontsize=14)
plt.ylim(0,1.05)
# plt.yticks( fontfamily='Sans serif', fontsize=14)
# plt.xlabel('Allele frequency', fontfamily='Sans serif', fontsize=14)
plt.ylabel('Frequency', fontfamily='Sans serif', fontsize=14)
plt.title(f'{timepoint}', fontfamily='Sans serif', fontsize=14)

for i, v in enumerate(normalized_counts):
    plt.text(i, v + 0.02, f'{v:.3f}({hets_counts[1][i]})', ha='center', fontfamily='Sans serif', fontsize=12)

plt.tight_layout()
plt.savefig(f'heterozygous{timepoint}.pdf', format='pdf', dpi=600)
plt.show()
plt.close()

