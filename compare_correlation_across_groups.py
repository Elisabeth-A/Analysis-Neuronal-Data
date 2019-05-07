#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  3 15:57:18 2019

@author: Elisabeth
"""
"""
Comparison of correlations across populations/groups (paired, unpaired, Striatum) using Fisher's Z statistic

"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import seaborn as sns
from scipy.stats import spearmanr



# load data from matlab file which contains response size data for the different stimuli & timepoints and groups
Path='/Users/Elisabeth/Documents/python_scripts/data_AUC/'
filename= 'learners_ER_f5_ERreliable_cut_hab.mat'
filename_un= 'learners_ER_unpaired_f5_ERreliable_cut_hab.mat'
filename_str= 'learners_ER_striatum_f5_ERreliable_cut_hab.mat'

mat_pair=scipy.io.loadmat(Path+filename)
mat_un=scipy.io.loadmat(Path+filename_un)
mat_str=scipy.io.loadmat(Path+filename_str)

# combining CS1 and CS2 in unpaired group

un_plus=np.array(mat_un['AUC_plus'])
un_minus=np.array(mat_un['AUC_minus'])
comb_CS_un=np.hstack((un_minus,un_plus))
comb_CS_dic= {'AUC_comb':comb_CS_un}
mat_un.update(comb_CS_dic)

(hab_r, ret_r)=correl_lin_reg_groups(mat_pair['AUC_minus'], mat_un['AUC_comb'], mat_str['AUC_minus'])


# plotting correlation of responses of all groups in one plot 
def correl_lin_reg_groups(cond_mat1,cond_mat2, cond_mat3):
    
    # define subplots
    fig,(ax1,ax2)=plt.subplots(figsize=(10,5), ncols=2, nrows=1)
    
    # calculate spearman r and plot comparison 
    hab_r=[spearmanr(np.transpose(cond_mat1[0:2]))[0],
           spearmanr(np.transpose(cond_mat2[0:2]))[0],
           spearmanr(np.transpose(cond_mat3[0:2]))[0]]
    
    ret_r=[spearmanr(np.transpose(cond_mat1[1:3]))[0],
           spearmanr(np.transpose(cond_mat2[1:3]))[0],
           spearmanr(np.transpose(cond_mat3[1:3]))[0]]
              
    all_=np.hstack((np.hstack((cond_mat1,cond_mat2)), cond_mat3))         
        ## use seaborn library regplot 
    sns.regplot(x=np.transpose(cond_mat1[1]), y=np.transpose(cond_mat1[0]),ax=ax1, ci=None, color="r")
    sns.regplot(x=np.transpose(cond_mat2[1]), y=np.transpose(cond_mat2[0]),ax=ax1, ci=None, color="b")
    sns.regplot(x=np.transpose(cond_mat3[1]), y=np.transpose(cond_mat3[0]),ax=ax1, ci=None, color="g")
    
        
    ax1.legend(loc='upper left', labels=['paired', 'pseudo', 'striatum'] )
    ax1.set(xlabel='response size early recall', ylabel='response size habituation')
    ax1.set_xlim([np.floor(all_.min())-1,np.ceil(all_.max())+1])
    ax1.set_ylim([np.floor(all_.min())-1,np.ceil(all_.max())+1])

        ## use seaborn library regplot 
    sns.regplot(x=np.transpose(cond_mat1[1]), y=np.transpose(cond_mat1[2]),ax=ax2, ci=None, color="r")
    sns.regplot(x=np.transpose(cond_mat2[1]), y=np.transpose(cond_mat2[2]),ax=ax2, ci=None, color="b")
    sns.regplot(x=np.transpose(cond_mat3[1]), y=np.transpose(cond_mat3[2]),ax=ax2, ci=None, color="g")
        
    ax2.legend(loc='upper left', labels=['paired', 'pseudo', 'striatum'] )
    ax2.set(xlabel='response size early recall', ylabel='response size 24 h recall')
    ax2.set_xlim([np.floor(all_.min())-1,np.ceil(all_.max())+1])
    ax2.set_ylim([np.floor(all_.min())-1,np.ceil(all_.max())+1])
    
        
    # plotting unity line and set axis range
    # find min and max of all data
    ax1.plot(np.arange(np.floor(all_.min())-1,np.ceil(all_.max())+1,1),
                np.arange(np.floor(all_.min())-1,np.ceil(all_.max())+1,1),'k--')
 
    ax2.plot(np.arange(np.floor(all_.min())-1,np.ceil(all_.max())+1,1),
             np.arange(np.floor(all_.min())-1,np.ceil(all_.max())+1,1),'k--')   
        
    return hab_r, ret_r



    # Fisher z statistic to compare correlation coefficient across groups
    # input are the correlation coeffiecients (r1,r2) and the number of samples (N1,N2)
def FisherZ_corr(r1, r2, N1, N2):    
    # transformation of correlation coefficients
   r_trans1=(0.5)*np.log((1+r1)/(1-r1))
   r_trans2=(0.5)*np.log((1+r2)/(1-r2))

    #compute test statistic
   z=(r_trans1-r_trans2)/np.sqrt((1/(N1-3))+(1/(N2-3)))
   Z_to_p=scipy.stats.norm.sf(abs(z))
   
   print("p-value:", Z_to_p)


p_paired_pseudo=FisherZ_corr(hab_r[0], hab_r[1], 154, 198)
p_paired_str=FisherZ_corr(hab_r[0], hab_r[2], 154, 200)
p_pseudo_str=FisherZ_corr(hab_r[1], hab_r[2], 198, 200)
