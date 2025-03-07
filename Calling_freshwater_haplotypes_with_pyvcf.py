import vcf 
import numpy as np
from sys import argv
from matplotlib import rcParams
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import colors

'''
This function takes 4 arguments: 1: a file containing all the names of the vcf files. 2: a tab-delimited file with the SNPs in 
344 loci with 6 columns: chromosome name, the representative SNP for each locus, number of SNPs within the locus, positions matching the genome coordinates of the snps, freshwater alleles at these snps, marine (oceanic) alleles at these snps.
3: the absolute path to the vcf files given in argument 1
4: minimum number of SNPs that should be in a locus in order for haplotypes to be called for the locus (we set this to 3 in the paper).
'''

vcf_list_file=argv[1]# this contains a list of individual genomes in a vcf format 
snp_list_file=argv[2] #'fw_haps_poolSNPs.CH_SC_LB.2500_50.50000.phy_specific_BF.res'
path_to_vcf=argv[3] 
snps_min=argv[4]

file = open(vcf_list_file, 'r') #arg 2
vcf_path_data=file.read()
vcf_path_data=vcf_path_data.split('\n')
file.close()
if vcf_path_data[-1]=='':
    del(vcf_path_data[-1])
    
vcf_paths=[]
vcf_dir=path_to_vcf #arg3
for i in range(len(vcf_path_data)):
    #print(vcf_dir+vcf_path_data[i])
    vcf_paths.append(vcf_dir+vcf_path_data[i])

snp_file = open(snp_list_file, 'r') #arg1 
snp_data=snp_file.read()
snp_data=snp_data.split('\n')
snp_file.close()
if snp_data[-1]=='':
    del(snp_data[-1])

min_snps=snps_min 

def GLtoPB(GLs):
    PB=(10**GLs)/np.sum(10**GLs)
    return PB

All_samps_in_pop=[]
Samp_names=[]
rep_snps=[]
for k in range(len(vcf_paths)):
    #print(vcf_paths[k])
    vcf_reader=vcf.Reader(filename=vcf_paths[k]) #pick one sample vcf 
    samp_name=vcf_paths[k].split('/')[-1].split('.')[0]
    hap_geno=[]
    for i in range(1,len(snp_data)): #pick one tempopeak 
        #print(i)
        k=snp_data[i].split('\t')
        chrom=k[0]
        rep_snp=k[1]
        positions=k[3] #all positions within a tempopeak 
        fw_alleles=k[4].split(',')
        mar_alleles=k[5].split(',')
        pos=positions.split(',')
        geno_probs_samp=np.zeros((3,len(pos)),dtype='float32')
        geno_rec=[]
        for j in range(len(pos)):
            site=pos[j] #site within tempopeak 
            fw_allele=fw_alleles[j] #fw_allele at site
            mar_allele=mar_alleles[j] #marine allele at site 
            position_record_vcf = vcf_reader.fetch(chrom, int(site)-1, int(site)) #extract the tempo peak snp from vcf  
            
            for record in position_record_vcf:
                f=record.samples
                l=record.POS
                gt=record.samples[0]['GT']
                ref=record.REF
                alts=record.ALT
                gp=record.samples[0]['GP']
                gen_vals=gt.split('|')

                if len(alts)==1 and fw_allele == ref and mar_allele == alts[0]: #only considering biallelic sites 
                    geno_rec.append(gp)
                elif len(alts)==1 and fw_allele == alts[0]  and mar_allele ==ref :
                    gp_sam =gp[2],gp[1],gp[0]
                    geno_rec.append(gp_sam)
        hap_geno.append(geno_rec) #record of haplotype of all tempo_peaks for one individual 344 total 
    rep_snps.append(rep_snp)
    All_samps_in_pop.append(hap_geno)
    Samp_names.append(samp_name)
    

ma_count=np.zeros((len(vcf_path_data),len(snp_data)-1),dtype='int32')
probs=np.zeros((len(vcf_path_data),len(snp_data)-1,3),dtype='float32')
ma_count[:]=-9
probs[:]=-9

for k in range(len(All_samps_in_pop)):
    tempo_samp_data=All_samps_in_pop[k]
    for i in range(len(tempo_samp_data)):
        if len(tempo_samp_data[i]) > min_snps:
            hap_l=np.sum(np.array(tempo_samp_data[i]),axis=0)
            PB=GLtoPB(np.clip(hap_l-np.max(hap_l),-300,0))
            ma_count[k][i]=np.argmax(PB)


mask=ma_count!=-9
ma_count=ma_count[mask].reshape(ma_count.shape[0], -1)
indices=np.where(ma_count[0] != -9)

k=snp_data[1].split('\t')
chromo=k[0]
chrom_lines=[0]
print(chrom_lines)
chrom_order=[]
chrom_lines_order=[0]
for g in range(1,len(indices[0])):
    k=snp_data[g].split('\t')
    if g in indices[0] and k[0]!=chromo:
        chrom_lines.append(g-1)
        chrom_order.append(chromo)
        chromo=k[0]
        chrom_lines_order.append([chromo, g-1])
        print(chromo, g)

chrom_lines.append(g)
    
chrom_order.append(chromo)
chrom_lines_order.append([chromo, g])  

chrom_lines_order[0]=['chrI', 0]

Rom_to_int = {'I': 1, 'II': 2, 'III': 3, 'IV': 4, 'V': 5, 'VI': 6, 'VII': 7, 'VIII': 8, 'IX': 9, 'X': 10, 'XI': 11, 'XII': 12, 'XIII': 13, 'XIV': 14, 'XIX': 19, 'XVI': 16, 'XVIII': 18, 'XX': 20, 'XXI': 21}

chr_sizes=[]
for j in range(1,len(chrom_lines)):
    chr_sizes.append(chrom_lines[j]- chrom_lines[j-1])
    
chrom_order_line_new=[]
for i in range(1,len(chrom_order)):
    chrom_order_line_new.append([chrom_lines_order[i-1][0], chr_sizes[i-1], chrom_lines_order[i-1][1], Rom_to_int[chrom_order[i-1][3:]]])
    
    
chrom_order_line_new.append([chrom_lines_order[i][0], chr_sizes[i],chrom_lines_order[i][1], Rom_to_int[chrom_order[i][3:]]])

chrom_order_line_new_sort=sorted(chrom_order_line_new, key=lambda x: x[3])  


new_chrom_order=[]
new_chrom_lines=[0]
chrom_order_line_new_sort
chrom=0
for a in range(len(chrom_order_line_new_sort)):
    chrom=chrom_order_line_new_sort[a][1]+chrom
    new_chrom_order.append(chrom_order_line_new_sort[a][0])
    new_chrom_lines.append(chrom)


fileout=filein+str(len(indices))    

pp = PdfPages(fileout+'.pdf')

rcParams['figure.figsize'] = 0.1*(len(indices[0])-1), len(vcf_path_data)*0.05



fig, ax = plt.subplots()

if np.max(ma_count_new)==1:
    if np.min(ma_count_new)==0:
        cmap_new = colors.ListedColormap(['#155084','#edea18']) #blue, yellow,
    elif np.min(ma_count_new)==-9:
        cmap_new = colors.ListedColormap(['#767676','#155084','#edea18']) #grey, blue, yellow,
    else:
        print('aaaa')
elif np.max(ma_count_new)==2:
    if np.min(ma_count_new)==0:
        cmap_new = colors.ListedColormap(['#155084','#edea18','#9f2305'])  #blue, yellow, red
    elif np.min(ma_count_new)==-9:
        cmap_new = colors.ListedColormap(['#767676','#155084','#edea18','#9f2305'])  #grey, blue, yellow, red
    else:
        print('aaaa')
cmap = ax.pcolormesh(np.clip(np.flipud(ma_count_new),-1,2),cmap=cmap_new,linewidth=0,rasterized=True)



ax.set_yticks([0,len(vcf_path_data)])
ax.set_yticklabels('')


ax.set_xticks(np.arange(1,len(indices[0]))-0.5,minor=False)
ax.set_xticklabels('')

ax.set_xlim(0, len(indices[0])-1)
ax.set_ylim(0, len(vcf_path_data))

for g in range(len(chrom_lines)):
    plt.axvline(x=chrom_lines[g],ymin=-0.025,ymax=1.025,clip_on=False,color="black")

for g in range(1,len(chrom_lines)):
    plt.text(((chrom_lines[g]-chrom_lines[g-1])/2.0)+chrom_lines[g-1], -2, chrom_order[g-1][3:], fontsize=8,ha='center',va='center')#, transform=plt.gcf().transFigure)

plt.savefig(pp, format='pdf', bbox_inches=0, dpi=600)
plt.close()

pp.close()

    



