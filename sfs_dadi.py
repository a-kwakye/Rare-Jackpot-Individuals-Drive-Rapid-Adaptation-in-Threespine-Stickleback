import dadi
import nlopt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def read_unfoldedsfs(sfs):
    f=40
    k=sfs.split('/')
    if len(k)> 1:
        pop=k[-1]
    else:
        pop=k
        
    data = np.genfromtxt(sfs, delimiter=' ')
    out_data=np.sum(data)
    spectrum=dadi.Spectrum(data, pop_ids=[pop])
    
    n=int((len(spectrum)-1)/2)
    pi=spectrum.pi()/np.sum(data)
    Wat_theta=spectrum.Watterson_theta()/np.sum(data)
    Taj=spectrum.Tajima_D()
    
    fold_spectrum=spectrum.fold()
    
    div=[pi, Wat_theta, Taj]
    proj=spectrum.project([f])
    pi_proj=proj.pi()/np.sum(data)
    Wat_theta_proj=proj.Watterson_theta()/np.sum(data)
    Taj_proj=proj.Tajima_D()
    div_proj=[pi_proj, Wat_theta_proj, Taj_proj]
    return fold_spectrum ,  div_proj, out_data


def num_of_Sites(sfs_dir, sfs_list,timepoint):
    file_sfs=sfs_dir+sfs_list
    file=open(file_sfs, 'r')
    data=file.read()
    data_out=data.split('\n')
    if data_out[-1]=="":
        del(data_out[-1])
    timpoint_sites=[]
    for i in data_out:
        k=i.split('_')  
        
        if k[0] == timepoint:
            
            a=read_unfoldedsfs(sfs_dir+i)[2]
            timpoint_sites.append(a)
    return timpoint_sites


def load_sfs(sfs_dir, sfs_list,timepoint):
    file_sfs=sfs_dir+sfs_list
    file=open(file_sfs, 'r')
    data=file.read()
    data_out=data.split('\n')
    if data_out[-1]=="":
        del(data_out[-1])
        
    timepoint_sfs=[]
    timepoint_div=[]
    
    for i in data_out:
        k=i.split('_')  
        
        if k[0] == timepoint:
            timepoint_sfs.append(read_unfoldedsfs(sfs_dir+i)[0])
            timepoint_div.append(read_unfoldedsfs(sfs_dir+i)[1])
    return timepoint_sfs, timepoint_div




n=40

sfs_dir='/gpfs/scratch/akwakye/AIM_1/transposons_masked_ref/GL_2/realsfs_sfs/unfolded/'
#unfolded_list='list_unfolded_sfs'
unfolded_list='list_unfolded_high_sfs'


RS2019_sfs=load_sfs(sfs_dir, unfolded_list, 'RS2019')


SC2013_sfs=load_sfs(sfs_dir, unfolded_list, 'SC2013')

SC2014_sfs=load_sfs(sfs_dir, unfolded_list, 'SC2014')
SC2015_sfs=load_sfs(sfs_dir, unfolded_list, 'SC2015')
SC2017_sfs=load_sfs(sfs_dir, unfolded_list, 'SC2017')

SC2020_sfs=load_sfs(sfs_dir, unfolded_list, 'Scout2020')
#SC2020_sfs=load_sfs(sfs_dir, unfolded_list, 'sc2020')


RS2019_sfs_projected=project_sfs(RS2019_sfs[0], n)

SC2013_sfs_projected=project_sfs(SC2013_sfs[0], n)

SC2014_sfs_projected=project_sfs(SC2014_sfs[0], n)

SC2015_sfs_projected=project_sfs(SC2015_sfs[0], n)

SC2017_sfs_projected=project_sfs(SC2017_sfs[0], n)



def barplot_foldedsfs(zero_fold, four_fold, WG, k, ax):
    fileout = zero_fold.pop_ids[0].split('_')[0]
    if fileout =='Scout2020':
        fileout='SC2020'
    print(fileout)
    to_plot_0 = zero_fold[1:k]
    to_plot_4 = four_fold[1:k]
    to_plot_wg = WG[1:k]

    bar_width = 0.30
    x = np.arange(0, k-1)

    # Plot bars on the provided axis (ax)
    ax.bar(x - bar_width, to_plot_wg / np.sum(WG), width=bar_width, color='#9F2305', label='Whole genome')
    #print(WG.pop_ids)
    ax.bar(x, to_plot_0[0:] / np.sum(zero_fold), width=bar_width, label='0-fold', color='#155084', alpha=0.7)
    #print(zero_fold.pop_ids)
    ax.bar(x + bar_width, to_plot_4[0:] / np.sum(four_fold), width=bar_width, label='4-fold', color='#819D4E', alpha=0.7)
    #print(four_fold.pop_ids)
    ax.set_xticks(range(0, k-1))
    ax.set_xticklabels(list(range(1, k)),  fontfamily='Sans serif', fontsize=7)
    
    ax.set_xlabel('Minor allele frequency', fontfamily='Sans serif', fontsize=14)
    ax.set_ylabel('Proportion', fontfamily='Sans serif', fontsize=14)
    ax.set_title(f' {fileout} ', fontfamily='Sans serif', fontsize=14)
    
    # Add legend
    ax.legend()

def create_combined_plot(k, RS2019_sfs, SC2013_sfs, SC2014_sfs, SC2015_sfs, SC2017_sfs):
    
    # Create a figure with 6 subplots arranged in a 3x2 grid
    fig, axs = plt.subplots(3, 2, figsize=(12, 12))  # Adjust figure size as needed
    
    # Flatten the axs array to make it easier to index
    axs = axs.flatten()

    # Labels for the panels
    panel_labels = ['A', 'B', 'C', 'D', 'E', 'F']

    # First, generate all subplots and store the maximum y-value across all subplots for setting uniform ylim
    max_y = 0  # Variable to track the maximum y-value across all subplots

    # Generate the subplots and add panel labels
    barplot_foldedsfs(RS2019_sfs[0][0], RS2019_sfs[0][1], RS2019_sfs[0][2], k, ax=axs[0])
    axs[0].text(-0.1, 1.05, 'A', transform=axs[0].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
    max_y = max(max_y, axs[0].get_ylim()[1])  # Update max_y

    barplot_foldedsfs(SC2013_sfs[0][0], SC2013_sfs[0][1], SC2013_sfs[0][2], k, ax=axs[1])
    axs[1].text(-0.1, 1.05, 'B', transform=axs[1].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
    max_y = max(max_y, axs[1].get_ylim()[1])

    barplot_foldedsfs(SC2014_sfs[0][0], SC2014_sfs[0][1], SC2014_sfs[0][2], k, ax=axs[2])
    axs[2].text(-0.1, 1.05, 'C', transform=axs[2].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
    max_y = max(max_y, axs[2].get_ylim()[1])

    barplot_foldedsfs(SC2015_sfs[0][0], SC2015_sfs[0][1], SC2015_sfs[0][2], k, ax=axs[3])
    axs[3].text(-0.1, 1.05, 'D', transform=axs[3].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
    max_y = max(max_y, axs[3].get_ylim()[1])

    barplot_foldedsfs(SC2017_sfs[0][0], SC2017_sfs[0][1], SC2017_sfs[0][2], k, ax=axs[4])
    axs[4].text(-0.1, 1.05, 'E', transform=axs[4].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
    max_y = max(max_y, axs[4].get_ylim()[1])

    fig.delaxes(axs[5]) 
    # Set the same ylim for all subplots
    for ax in axs:
        ax.set_ylim(0, max_y)

    # Adjust layout to avoid overlap and make room for labels
    plt.tight_layout()

    # Save the figure
    if len(RS2019_sfs[0][0]) == len(SC2013_sfs[0][0]):
    
        plt.savefig(f"All_timepoints_WG_4_0_fold_dadi_folded_{k}_projected.pdf", format='pdf', dpi=600)
    elif len(RS2019_sfs[0][0]) != len(SC2013_sfs[0][0]):
        plt.savefig(f"All_timepoints_WG_4_0_fold_dadi_folded_{k}.pdf", format='pdf', dpi=600)

    # Show the figure
    plt.show()


def barplot_sfs_all_WG(RS2019_sfs, SC2013_sfs, SC2014_sfs, SC2015_sfs, SC2017_sfs, SC2020_sfs, k, bar_width, n):
    to_plot_RS2019 = RS2019_sfs[1:k]
    to_plot_SC2013 = SC2013_sfs[1:k]
    to_plot_SC2014 = SC2014_sfs[1:k]
    to_plot_SC2015 = SC2015_sfs[1:k]
    to_plot_SC2017 = SC2017_sfs[1:k]
    to_plot_SC2020 = SC2020_sfs[1:k]
    
    x = np.arange(0, k-1)
    textsize=18
    # Set figure size (width 336mm, height adjusted to maintain proportions)
    plt.figure(figsize=(336/25.4, 210/25.4))  
    
    plt.bar(x - bar_width*2, to_plot_RS2019, width=bar_width, color='#550000', label='RS2019') ##9F2305
    print(RS2019_sfs.pop_ids)
    
    plt.bar(x - bar_width, to_plot_SC2013, width=bar_width, color='#9F2305', label='SC2013') ##A93C08
    print(SC2013_sfs.pop_ids)
    
    plt.bar(x, to_plot_SC2014, width=bar_width, color='#EDEA18', label='SC2014')
    print(SC2014_sfs.pop_ids)
    
    plt.bar(x + bar_width, to_plot_SC2015, width=bar_width, label='SC2015', color='#819D4E')
    print(SC2015_sfs.pop_ids)
    
    plt.bar(x + bar_width*2, to_plot_SC2017, width=bar_width, label='SC2017', color='#155084')
    print(SC2017_sfs.pop_ids)

    plt.xticks(range(0, k-1), list(range(1, k)), rotation=90, fontfamily='Sans serif', fontsize=textsize)
    plt.yticks(fontfamily='Sans serif', fontsize=textsize)
    plt.xlabel('Minor allele frequency', fontfamily='Sans serif', fontsize=textsize)
    plt.ylabel('Count', fontfamily='Sans serif', fontsize=textsize)
    plt.title('Site frequency spectrum (whole genome) with angsd and dadi', fontfamily='Sans serif', fontsize=textsize)
    plt.tight_layout()
    
    plt.legend()
    
    if len(RS2019_sfs) == len(SC2014_sfs):
        plt.savefig(f'SFS_WG_dadi_folded_{k}_projected_count.pdf', format='pdf', dpi=600)
    else:
        plt.savefig(f'SFS_WG_dadi_folded_{n}_{k}_unprojected_count.pdf', format='pdf', dpi=600)

    plt.show()
    plt.close()

    

    
barplot_sfs_all_WG(RS2019_sfs_projected[0][2], SC2013_sfs_projected[0][2],SC2014_sfs_projected[0][2], SC2015_sfs_projected[0][2], SC2017_sfs_projected[0][2],SC2020_sfs[0][2], 21, 0.18, 21)


create_combined_plot(21, RS2019_sfs_projected, SC2013_sfs_projected, SC2014_sfs_projected, SC2015_sfs_projected,SC2017_sfs_projected)



barplot_sfs_all_WG(RS2019_sfs_projected[0][2], SC2013_sfs_projected[0][2],SC2014_sfs_projected[0][2], SC2015_sfs_projected[0][2], 
SC2017_sfs_projected[0][2],SC2020_sfs[0][2], 21, 0.18, 21)





