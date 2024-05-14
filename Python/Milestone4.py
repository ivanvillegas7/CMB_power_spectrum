# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 12:06:38 2024

@author: Iván Villegas Pérez
"""

#Import all relevant packages

import scipy as sc

import numpy as np

import matplotlib.pyplot as plt

import auxiliar as aux

def CMB_PowerSpectrum():
    
    data = np.loadtxt('../Results/cells.txt')
    
    ell: np.array(float) = data[:, 0]
    
    C_ell: np.array(float) = data[:, 1]
       
    plt.figure()
    plt.plot(ell, C_ell, label='Theory prediction')
    plt.xlabel(r'Multipole $\ell$')
    plt.ylabel(r'$C_\ell$')
    plt.title('CMB power-spectrum')
    plt.xscale('log')
    plt.legend()
    plt.grid()
    #plt.savefig('../Plots/Milestone IV/CMB_PS.pdf')
    
def Matter_PowerSpectrum(polarization: bool):
    
    data = np.loadtxt('../Results/cells.txt')
    
    ell: np.array(float) = data[:, 0]
    
    C_ell: np.array(float) = data[:, 1]
    
    data_gal = np.loadtxt('../Data/galaxy_survey_data.txt', skiprows=1)
    
    k_gal: np.array(float) = data_gal[:, 0]
    
    P_gal: np.array(float) = data_gal[:, 1]
    
    ErrorP_gal: np.array(float) = data_gal[:, 2]
    
    data_ACT = np.loadtxt('../Data/WMAP_ACT_data.txt', skiprows=1)
    
    k_ACT: np.array(float) = data_ACT[:, 0]
    
    P_ACT: np.array(float) = data_ACT[:, 1]
    
    P_upper: np.array(float) = data_ACT[:, 2]
    
    data_lowTT = np.loadtxt('../Data/low_TT.txt', skiprows=1)
    
    l_lowTT: np.array(float) = data_lowTT[:, 0]
    
    D_l_lowTT: np.array(float) = data_lowTT[:, 1]
    
    DeltaD_down_lowTT: np.array(float) = data_lowTT[:, 2]
    
    DeltaD_up_lowTT: np.array(float) = data_lowTT[:, 3]
    
    if polarization:
        
        data_highTT = np.loadtxt('../Data/high_TT.txt', skiprows=1)
        
        l_highTT: np.array(float) = data_highTT[:, 0]
        
        D_l_highTT: np.array(float) = data_highTT[:, 1]
        
        DeltaD_down_highTT: np.array(float) = data_highTT[:, 2]
        
        DeltaD_up_highTT: np.array(float) = data_highTT[:, 3]
    
        data_highEE = np.loadtxt('../Data/high_EE.txt', skiprows=1)
        
        l_highEE: np.array(float) = data_highEE[:, 0]
        
        D_l_highEE: np.array(float) = data_highEE[:, 1]
        
        DeltaD_down_highEE: np.array(float) = data_highEE[:, 2]
        
        DeltaD_up_highEE: np.array(float) = data_highEE[:, 3]
        
        data_highTE = np.loadtxt('../Data/high_TE.txt', skiprows=1)
        
        l_highTE: np.array(float) = data_highTE[:, 0]
        
        D_l_highTE: np.array(float) = data_highTE[:, 1]
        
        DeltaD_down_highTE: np.array(float) = data_highTE[:, 2]
        
        DeltaD_up_highTE: np.array(float) = data_highTE[:, 3]
        
    plt.figure()
    #plt.plot(ell, ell*(ell+1)*C_ell/(2*np.pi), label='Theory prediction')
    plt.errorbar(l_lowTT, D_l_lowTT, yerr=[DeltaD_down_lowTT, DeltaD_up_lowTT],\
                 ls='none', label=r'Low $\ell$ TT data', marker='.', capsize=2)
    if polarization:
        plt.errorbar(l_highTT, D_l_highTT,\
                     yerr=[DeltaD_down_highTT, DeltaD_up_highTT], ls='none',\
                     label=r'High $\ell$ TT data', marker='.', capsize=2)
        plt.errorbar(l_highEE, D_l_highEE,\
                     yerr=[DeltaD_down_highEE, DeltaD_up_highEE], ls='none',\
                     label=r'High $\ell$ EE data', marker='.', capsize=2)
        plt.errorbar(l_highTE, D_l_highTE,\
                     yerr=[DeltaD_down_highTE, DeltaD_up_highTE], ls='none',\
                     label=r'High $\ell$ TE data', marker='.', capsize=2)
    plt.xlabel(r'Multipole $\ell$')
    plt.ylabel(r'$\ell(\ell+1)C_\ell/2\pi$ [$\mu$K$^2$]')
    plt.title('Matter power-spectrum')
    plt.xscale('log')
    #plt.yscale('log')
    plt.legend()
    plt.grid()
    #plt.savefig('../Plots/Milestone IV/Matter PS.pdf')
    
    #Read the data from 'cosmology.txt' in folder 'Results'.
    
    data = np.loadtxt('../Results/cosmology.txt')
    
    #Hp: conformal Hubble factor.
    
    Hp: np.array(float) = data[:, 2]
    
    #Get index of radiation and matter equality.
    
    index: int = aux.index_equality()[0]
    
    Hp_eq: float = (Hp[index]/(100*1e3/(1e6*sc.constants.parsec)))
    
    k_eq: float = Hp_eq*1e5/(sc.constants.c*sc.constants.h)
        
    plt.figure()
    plt.errorbar(k_gal, P_gal, ErrorP_gal, label='SDSS Galaxies (DR7 LRG)',\
                 ls='none', marker='.', capsize=2)
    plt.errorbar(k_ACT, P_ACT, P_upper, label='CMB (WMAP+ACT)', ls='none',\
                 marker='.', capsize=2)
    plt.vlines(k_eq, plt.ylim()[0], plt.ylim()[1], label=r'$k_\text{eq}$',\
               ls='--')
    plt.xlabel(r'Wavenumber $k$ [$h$/Mpc]')
    plt.ylabel(r'$P(k)$ [(Mpc/$h$)$^2$]')
    plt.title('The total matter power-spectrum')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.grid()
    #plt.savefig('../Plots/Milestone IV/Total matter PS.pdf')
    
    
def milestone4(polarization: bool):
    
    CMB_PowerSpectrum()
    Matter_PowerSpectrum(polarization)
    
milestone4(False)
    