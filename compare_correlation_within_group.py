#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 11:30:33 2019

@author: Elisabeth
"""

"""
correlation analysis on neuronal response size data 
comparing different timepoints (recording sessions: habituation, early recall and 24 h recall)
and for each stimulus (plus and minus)

"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import seaborn as sns
from scipy.stats import spearmanr

# load data from matlab file which contains response size data for the different stimuli & 3 timepoints
Path='/Users/Elisabeth/Dropbox/python_scripts/data_AUC/'
filename= 'learners_ER_f5_ERreliable_cut_hab.mat'


mat_pair=scipy.io.loadmat(Path+filename)


## calculate correlation coefficient for each stimulus (minus or plus)

def main():
    paired_r=correl_lin_reg(mat_pair['AUC_minus'])
    
    Steiger_z(paired_r)

    p_KS_minus=bootstrap_r(paired_r,mat_pair['AUC_minus'])  

    return paired_r,p_KS_minus

def correl_lin_reg(cond_mat):
    # calculate spearman r for 3 different time points (column in input matrix)
    paired_r=[spearmanr(np.transpose(cond_mat[0:2]))[0],
                  spearmanr(np.transpose(cond_mat[1:3]))[0],
                  spearmanr(np.transpose(cond_mat[[0,2],:]))[0]]
              


    ## use seaborn library regplot to plot data with regression line
    f, ax = plt.subplots(figsize=(5, 5))
    sns.regplot(x=np.transpose(cond_mat[1]), y=np.transpose(cond_mat[0]))
    sns.regplot(x=np.transpose(cond_mat[1]), y=np.transpose(cond_mat[2]))
    
    plt.legend(loc='upper left', labels=['habituation', '24h recall'] )
    plt.ylabel('response size')
    plt.xlabel('response size early recall')
    plt.axis([np.floor(cond_mat.min()),np.ceil(cond_mat.max()), np.floor(cond_mat.min()),np.ceil(cond_mat.max())])
    # plotting unity line
    plt.plot(np.arange(np.floor(cond_mat.min()),np.ceil(cond_mat.max())+1,1),
            np.arange(np.floor(cond_mat.min()),np.ceil(cond_mat.max())+1,1),'k--')
    
    
    return paired_r


    
def Steiger_z(paired_r):
    # compare correlation coefficient across time points within each group (correlated samples)
    # use Steiger's z-test for "correlated correlations" within a population
    rm_sq=(paired_r[0]**2+paired_r[2]**2)/2
    f=(1-paired_r[2])/(2*(1-rm_sq))
    h=(1-(f*rm_sq))/(1-rm_sq)
    # calculate z transform of paired_r using arctanh
    paired_r_z=np.arctanh(paired_r)
    # Z-score and corresponding p-value
    Z = (paired_r_z[0] -paired_r_z[1])*(np.sqrt( len(mat_pair['AUC_minus'][0])-3)/np.sqrt(2*(1-paired_r[2])*h))
    Z_to_p=scipy.stats.norm.sf(abs(Z)) # one-tailed
    
    print('Steigers Z, p-value:',Z, Z_to_p)
    
def bootstrap_r(paired_r, cond_mat):
    # use bootstrapping to compare correlations (5000 iterations)
    boot_r_ER_24h=np.array([])
    boot_r_ER_hab=np.array([])
    
    for it in range(5000):
        subsample = np.random.choice(np.arange(0,len(cond_mat[1]),1),size = len(cond_mat[1]), replace = True)
        boot_r_ER_24h=np.append(boot_r_ER_24h, spearmanr(np.transpose(cond_mat[1:3])[subsample])[0])
        boot_r_ER_hab=np.append(boot_r_ER_hab, spearmanr(np.transpose(cond_mat[0:2])[subsample])[0])
    
    f, ax = plt.subplots(figsize=(5, 5))
    plt.hist(boot_r_ER_24h, bins='auto',align="left")
    plt.axvline(paired_r[1], color='b')
    
    plt.hist(boot_r_ER_hab, bins='auto',align="left")
    plt.axvline(paired_r[0], color='y')
    plt.xlabel('r')
    plt.ylabel('N occurences')
   
    #plt.xlim([np.floor(boot_r_ER_24h.min())-1,np.ceil(boot_r_ER_24h.max())+1])    
        
    #perform Kolomogorov Smirnoff test
    p_KS=scipy.stats.ks_2samp(boot_r_ER_24h, boot_r_ER_hab)
    
    print('KS p-value',p_KS[1])
  

if __name__ == '__main__':
    # in that case - call the main function
    main()
