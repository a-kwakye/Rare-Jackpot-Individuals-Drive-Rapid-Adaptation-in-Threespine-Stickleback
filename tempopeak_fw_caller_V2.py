from sys import argv
import pysam
import math
import numpy as np
import string
import numpy.ma as ma
import random
import matplotlib
matplotlib.use('Agg')
from matplotlib import rcParams
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import colors
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.optimize import fmin
from scipy.optimize import fminbound
from scipy.optimize import curve_fit
import time
from random import randint
from random import shuffle
import copy

bam_list_file=argv[1]#'RS2019_bam.list'#
snp_list_file=argv[2] #'fw_haps_poolSNPs.CH_SC_LB.2500_50.50000.phy_specific_BF.res'#argv[2]
ref_file='/vault/veeramah/people/kwakye/GasAcu1-4_genome/gasAcu1-4.fa'#argv[3]
filenameout=bam_list_file+'_fw_poolSNPs'
min_RD=1
MQ=30
BQ=20

min_snps=3 #minimum number of snps with read data required to call a haplotyppe

all_dic={}
all_dic['A']=0
all_dic['C']=1
all_dic['G']=2
all_dic['T']=3


###extract reads for a give position in a bam
def extract_bam_SNP(samfile,chromo,pos,BQ,MQ):
    var_list=[]
    for pileupcolumn in samfile.pileup(chromo,pos-1,pos,truncate=True,stepper='all'):
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                if (pileupread.alignment.mapping_quality>=MQ) and (ord(pileupread.alignment.qual[pileupread.query_position])-33>=BQ):
                        var_list.append([pileupread.alignment.query_sequence[pileupread.query_position],ord(pileupread.alignment.qual[pileupread.query_position])-33])
    return var_list


def phred2prob(x):
    return 10.0**(-x/10.0)

def prob2phred(x):
    return -10*math.log10(x)

###Phred_scale
def Phred_scale(ll):
    PL=-10*np.log10(ll/np.max(ll))
    PL=np.round(PL)
    PL=PL.astype(int)
    return PL

def GLtoPB(GLs):
    PB=(10**GLs)/np.sum(10**GLs)
    return PB

def geno_caller_3GT(X,ref,alt,all_dic):
    #diploid caller assuming that only assesses likelihood for three possible genotypes (ref/ref,ref/alt,alt/alt)
    GL=[0.0,0.0,0.0]

    count=0
    for g in range(len(X)):
        if all_dic.has_key(X[g][0])==False:
            continue
        err=phred2prob(X[g][1])*(1.0/3.0)
        tru=1-phred2prob(X[g][1])

        if X[g][0]==ref:
            GL[0]=GL[0]+math.log10(tru)
            GL[1]=GL[1]+math.log10((tru+err)/2)
            GL[2]=GL[2]+math.log10(err)
        elif X[g][0]==alt:
            GL[0]=GL[0]+math.log10(err)
            GL[1]=GL[1]+math.log10((tru+err)/2)
            GL[2]=GL[2]+math.log10(tru)
        else:
            GL[0]=GL[0]+math.log10(err)
            GL[1]=GL[1]+math.log10(err)
            GL[2]=GL[2]+math.log10(err)
        count+=1

    if count==0:
        GL=[0.0,0.0,0.0]
    return GL

###open up reference file
ref=pysam.FastaFile(ref_file)

###get bam.list
file = open(bam_list_file)
bam_list=file.read()
bam_list=string.split(bam_list,'\n')
file.close()

if bam_list[-1]=='':
    del(bam_list[-1])


###make a list of SNPs to be interogated
file = open(snp_list_file)
data=file.read()
data=string.split(data,'\n')
file.close()

if data[-1]=='':
    del(data[-1])

###open up bam files in dictionary (must be indexed)

bam_dic={}
samp_names=[]
for g in range(len(bam_list)):
    bam_dic[g] = pysam.AlignmentFile(bam_list[g], "rb")
    samp_names.append(string.split(bam_dic[g].header['RG'][0]['SM'],'/')[-1])


ma_count=np.zeros((len(samp_names),len(data)-1),dtype='int32')
probs=np.zeros((len(samp_names),len(data)-1,3),dtype='float32')
ma_count[:]=-9
probs[:]=-9


for g in range(1,len(data)):
    k=string.split(data[g],'\t')
    print str(g)+'_'+k[0]+'_'+k[1]
    chromo=k[0]
    positions=np.array(string.split(k[3],','),dtype='int32')
    fw_alleles=string.split(k[4],',')
    ma_alleles=string.split(k[5],',')

    for gg in range(len(bam_list)):
        genos=np.zeros((3,len(positions)),dtype='float32') #fw_homo,het,ma_homo
        snp_count=0
        for ggg in range(len(positions)):
            var_list=extract_bam_SNP(bam_dic[gg],chromo,positions[ggg],BQ,MQ)
            GL=geno_caller_3GT(var_list,fw_alleles[ggg],ma_alleles[ggg],all_dic)
            genos[:,ggg]=GL
            if len(var_list)>0:
                snp_count+=1

        if snp_count>=min_snps:
            hap_GL=np.sum(genos,axis=1)
            PB=GLtoPB(np.clip(hap_GL-np.max(hap_GL),-300,0))#GLtoPB(hap_GL)

            ma_count[gg][g-1]=np.argmax(PB)
            probs[gg][g-1]=PB

fileout=open(filenameout,'w')

for  g in range(len(ma_count)):
    ma_count_ind_nomiss=ma_count[g][np.where(ma_count[g]<>-9)[0]]
    fw_per_ind=np.sum(ma_count_ind_nomiss)/float(len(ma_count_ind_nomiss)*2)
    out=samp_names[g]+'\t'+str(fw_per_ind)+'\t'+str(len(ma_count_ind_nomiss))+'\n'
    fileout.write(out)

fileout.close()

k=string.split(data[1],'\t')
chromo=k[0]

chrom_lines=[0]
chrom_order=[]

for g in range(1,len(data)):
    k=string.split(data[g],'\t')
    if k[0]<>chromo:
        chrom_lines.append(g-1)
        chrom_order.append(chromo)
        chromo=k[0]

chrom_lines.append(g)
chrom_order.append(chromo)



pp = PdfPages(filenameout+'.pdf')

rcParams['figure.figsize'] = 0.1*(len(data)-1), len(samp_names)*0.05

fig, ax = plt.subplots()
#plt.title('Rabbit Slough '+chromo)
if np.max(ma_count)==1:
    if np.min(ma_count)==0:
        cmap_new = colors.ListedColormap(['#155084','#edea18']) #blue, yellow,
    elif np.min(ma_count)==-9:
        cmap_new = colors.ListedColormap(['#767676','#155084','#edea18']) #grey, blue, yellow,
    else:
        aaaa
elif np.max(ma_count)==2:
    #cmap_new = colors.ListedColormap(['#767676','#9f2305','#edea18','#155084'])
    if np.min(ma_count)==0:
        cmap_new = colors.ListedColormap(['#155084','#edea18','#9f2305'])  #blue, yellow, red
    elif np.min(ma_count)==-9:
        cmap_new = colors.ListedColormap(['#767676','#155084','#edea18','#9f2305'])  #grey, blue, yellow, red
    else:
        kjkj
cmap = ax.pcolormesh(np.clip(np.flipud(ma_count),-1,2),cmap=cmap_new,linewidth=0,rasterized=True)


ax.set_yticks([0,len(samp_names)])
ax.set_yticklabels('')


ax.set_xticks(np.arange(1,len(data))-0.5,minor=False)
ax.set_xticklabels('')

ax.set_xlim(0, len(data)-1)
ax.set_ylim(0, len(samp_names))

for g in range(len(chrom_lines)):
    plt.axvline(x=chrom_lines[g],ymin=-0.025,ymax=1.025,clip_on=False,color="black")

for g in range(1,len(chrom_lines)):
    plt.text(((chrom_lines[g]-chrom_lines[g-1])/2.0)+chrom_lines[g-1], -2, chrom_order[g-1][3:], fontsize=8,ha='center',va='center')#, transform=plt.gcf().transFigure)

plt.savefig(pp, format='pdf', bbox_inches=0, dpi=600)
plt.close()

pp.close()



#np.save('/gpfs/scratch/akwakye/AIM_1/FWALs_per_timepoint/SC_20_Samp_names.npy', Samp_names)
np.save('/gpfs/scratch/akwakye/AIM_1/FWALs_per_timepoint/SC_20_ma_count.npy', ma_count)
