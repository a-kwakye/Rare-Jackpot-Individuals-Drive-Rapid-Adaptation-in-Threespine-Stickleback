
import  Bio
from Bio import Phylo, AlignIO
import numpy as np
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import ClustalwCommandline
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib.pyplot as plt
import networkx as nx
import multiprocessing as mp
import os
from subprocess import Popen,PIPE
from sys import argv
import csv
import pysam


mito_dir='TN5_analyses/Mitochondrial_seq_TN5/'

List_of_files=open(mito_dir+'SC-14_fasta.list','r')
File_list=List_of_files.read()
SC_2014_Files=File_list.split('\n')
if SC_2014_Files[-1] =='':
    del(SC_2014_Files[-1])
    
List_of_files=open(mito_dir+'RS2019_fasta.list','r')
File_list=List_of_files.read()
RS2019_Files=File_list.split('\n')
if RS2019_Files[-1] =='':
    del(RS2019_Files[-1])
    
    
List_of_files=open(mito_dir+'SC-2013_fasta.list','r')
File_list=List_of_files.read()
SC_2013_Files=File_list.split('\n')
if SC_2013_Files[-1] =='':
    del(SC_2013_Files[-1])
       

List_of_files=open(mito_dir+'SC-2015_fasta.list','r')
File_list=List_of_files.read()
SC_2015_Files=File_list.split('\n')
if SC_2015_Files[-1] =='':
    del(SC_2015_Files[-1])


List_of_files=open(mito_dir+'SC-17_fasta.list','r')
File_list=List_of_files.read()
SC_2017_Files=File_list.split('\n')
if SC_2017_Files[-1] =='':
    del(SC_2017_Files[-1])

List_of_files=open(mito_dir+'SC2020.list','r')
File_list=List_of_files.read()
SC_2020_Files=File_list.split('\n')
if SC_2020_Files[-1] =='':
    del(SC_2020_Files[-1])

    


Files= RS2019_Files+ SC_2013_Files+SC_2014_Files+SC_2015_Files+SC_2017_Files+SC_2020_Files
 




ref_file='ChrM.fasta'
mito_ref=SeqIO.read(ref_file, "fasta")
mito_ref_seq=mito_ref.seq

seq_array=np.full((len(Files),16500),None, dtype=object)
#compare mito sequences with reference only output non-similar sites 
Seg_sites=[]
Seg_sites_samps=[]
samp_names=[]
nucleotides=["A","C", "G","T"]
for i in range(len(Files)):
    Sample_path= Files[i]
    samp_seg_sites=[]
    samp_sites=[]
    for record in SeqIO.parse(Sample_path, "fasta"): #takes one mitochondrial sequence 
        Mito=record.seq[0:16500]
        samp_names.append(record.id)
        #print(record.id)
        #print(len(mito_ref_seq), len(Mito))
        seq_array[i]=Mito

seq_array_t=np.transpose(seq_array)
samp_names_ar=np.array(samp_names)
min_count=22.5
min_freq=min_count/len(Files)

max_freq=1-min_freq
sites_filt=[]
nucleotides=["A","C", "G","T", 'N']
non_conform=[]
for i in range(len(seq_array_t)):
    k=np.unique(seq_array_t[i])
    kk=np.unique(seq_array_t[i], return_counts=True)
    if len(k) > 1 and len(k) < 3:
        if all(str.upper(e) in nucleotides for e in k):
            if str.upper(k[0]) != str.upper(k[1]) or str.upper(k[1]) != str.upper(k[0]):
                g=kk[1]/len(Files)
                
                if g[0] > g[1]:
                    if g[0] < max_freq:
                        #print(g)
                        #print(g,i)
                        sites_filt.append(i)
                elif g[0] < g[1]:
                    if g[0] > min_freq:
                        #print(g)
                        sites_filt.append(i)
                
#     elif len(k) > 2:
#         if any(str.upper(e) not in nucleotides for e in kk[0]):
#             indices_non=[iii for iii, e in enumerate(kk[0]) if str.upper(e) not in nucleotides]
#             indices_act=[iii for iii, e in enumerate(kk[0]) if str.upper(e) in nucleotides]
        
#             if len(kk[0]) - len(indices_non) ==2:
#                 g=kk[1][indices_act]/sum(kk[1][indices_act])
#                 #print(indices_act, indices_non)
#                # print(g, kk[0][indices_act], kk[1][indices_act])
#                 if g[0] > g[1]:
#                     if g[0] < max_freq:
#                         sites_filt.append(i)
#                 elif g[0] < g[1]:
#                     if g[0] > min_freq:
#                         sites_filt.append(i)   


