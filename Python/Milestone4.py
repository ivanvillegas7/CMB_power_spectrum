# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 12:06:38 2024

@author: Iván Villegas Pérez
"""

#Import all relevant packages

import scipy as sc

import numpy as np

import matplotlib.pyplot as plt

import healpy as hp

import auxiliar as aux

def CMB_PowerSpectrum():
    
    """
    Plot the Cosmic Microwave Background (CMB) power-spectrum.

    This function reads data from a file containing CMB power-spectrum
    information, including multipole moments (ell) and their corresponding
    power values (C_ell). It then generates a plot to visualize the theoretical
    prediction of the CMB power-spectrum.
    
    Parameters:
        None.

    Returns:
        None.
    """
    
    #Read the data from 'cells.txt' in folder 'Results'.
    
    data = np.loadtxt('../Results/cells.txt')
    
    #ell: multipole moment.
    
    ell: np.array(float) = data[:, 0]
    
    #C_ell_TT: power-spectrum with only the temperature source function.
    
    C_ell_TT: np.array(float) = data[:, 1]
    
    #Make the plot.
       
    plt.figure()
    plt.plot(ell, C_ell_TT, label='Theory prediction')
    plt.xlabel(r'Multipole $\ell$')
    plt.ylabel(r'$\frac{\ell(\ell+1)}{2\pi}C_\ell$ [$\mu$K$^2$]')
    plt.title('CMB power-spectrum')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.grid()
    plt.savefig('../Plots/Milestone IV/CMB_PS.pdf')
    
def Matter_PowerSpectrum(polarization: bool):
    
    """
    Plot the matter power spectrum.

    This function reads data from files containing matter power-spectrum
    information, including multipole moments (ell) and their corresponding
    power values (C_ell), as well as data from various surveys and experiments.
    It then generates a plot to visualize the theoretical prediction of the
    matter power-spectrum along with observational data.

    Parameters:
        polarization (bool): indicates if polarization has been included.

    Returns:
        None.
    """
    
    #Read the data from 'cells.txt' and 'Matter_PS.txt' in folder 'Results'
    
    data_Cells = np.loadtxt('../Results/cells.txt')
    
    #ell: multipole moment.
    
    ell: np.array(float) = data_Cells[:, 0]
    
    #C_ell: power-spectrum.
    
    C_ell_TT: np.array(float) = data_Cells[:, 1]
    
    data_MPS = np.loadtxt('../Results/Matter_PS.txt')
    
    #k: wavenumber.
    
    k: np.array(float) = data_MPS[:, 0]
    
    #Matter_PS: matter power-spectrum.
    
    Matter_PS: np.array(float) = data_MPS[:, 1]
    
    #Read and assign the data from 'galaxy_survey_data.txt' in folder 'Data'.
    
    data_gal = np.loadtxt('../Data/galaxy_survey_data.txt', skiprows=1)
    
    k_gal: np.array(float) = data_gal[:, 0]
    
    P_gal: np.array(float) = data_gal[:, 1]
    
    ErrorP_gal: np.array(float) = data_gal[:, 2]
    
    #Read and assign the data from 'WMAP_ACT_data.txt' in folder 'Data'.
    
    data_ACT = np.loadtxt('../Data/WMAP_ACT_data.txt', skiprows=1)
    
    k_ACT: np.array(float) = data_ACT[:, 0]
    
    P_ACT: np.array(float) = data_ACT[:, 1]
    
    P_upper: np.array(float) = data_ACT[:, 2]
    
    #Read and assign the data from 'low_TT.txt' in folder 'Data'.
    
    data_lowTT = np.loadtxt('../Data/low_TT.txt', skiprows=1)
    
    l_lowTT: np.array(float) = data_lowTT[:, 0]
    
    D_l_lowTT: np.array(float) = data_lowTT[:, 1]
    
    DeltaD_down_lowTT: np.array(float) = data_lowTT[:, 2]
    
    DeltaD_up_lowTT: np.array(float) = data_lowTT[:, 3]
    
    #Read and assign the data from 'high_TT.txt' in folder 'Data'.
    
    data_highTT = np.loadtxt('../Data/high_TT.txt', skiprows=1)
    
    l_highTT: np.array(float) = data_highTT[:, 0]
    
    D_l_highTT: np.array(float) = data_highTT[:, 1]
    
    DeltaD_down_highTT: np.array(float) = data_highTT[:, 2]
    
    DeltaD_up_highTT: np.array(float) = data_highTT[:, 3]
    
    if polarization:
        
        #Read and assign the data from 'high_EE.txt' in folder 'Data'.
    
        data_highEE = np.loadtxt('../Data/high_EE.txt', skiprows=1)
        
        l_highEE: np.array(float) = data_highEE[:, 0]
        
        D_l_highEE: np.array(float) = data_highEE[:, 1]
        
        DeltaD_down_highEE: np.array(float) = data_highEE[:, 2]
        
        DeltaD_up_highEE: np.array(float) = data_highEE[:, 3]
        
        #Read and assign the data from 'high_TE.txt' in folder 'Data'.
        
        data_highTE = np.loadtxt('../Data/high_TE.txt', skiprows=1)
        
        l_highTE: np.array(float) = data_highTE[:, 0]
        
        D_l_highTE: np.array(float) = data_highTE[:, 1]
        
        DeltaD_down_highTE: np.array(float) = data_highTE[:, 2]
        
        DeltaD_up_highTE: np.array(float) = data_highTE[:, 3]
        
    #Read and assign the data from 'Ly_alpha.txt' in folder 'Data'.
    
    data_Ly = np.loadtxt('../Data/Ly_alpha.txt', skiprows=1)
    
    k_Ly: np.array(float) = data_Ly[:, 0]
    
    P_Ly: np.array(float) = data_Ly[:, 1]
    
    DeltaP_Ly: np.array(float) = data_Ly[:, 2]
        
    #Make the first plot.
        
    plt.figure()
    plt.plot(ell, C_ell_TT, label='Theory prediction')
    plt.errorbar(l_lowTT, D_l_lowTT, yerr=[DeltaD_down_lowTT, DeltaD_up_lowTT],\
                 ls='none', label=r'Low $\ell$ TT data', marker='.', capsize=2)
    plt.errorbar(l_highTT, D_l_highTT,\
                 yerr=[DeltaD_down_highTT, DeltaD_up_highTT], ls='none',\
                 label=r'High $\ell$ TT data', marker='.', capsize=2)
    if polarization:
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
    plt.yscale('log')
    plt.legend()
    plt.grid()
    plt.savefig('../Plots/Milestone IV/Matter PS.pdf')
    
    #Read the data from 'cosmology.txt' in folder 'Results'.
    
    data = np.loadtxt('../Results/cosmology.txt')
    
    #Get index of radiation and matter equality.
    
    index: int = aux.index_equality()[0]
    
    #Hp_eq: conformal Hubble factor when there is radiation and matter equality.
    
    Hp_eq: np.array(float) = data[index, 2]*10*sc.constants.parsec #100km/(Mpc*s)
    
    #k_eq: wavenumber when there is radiation and matter equality.
    
    k_eq: float = Hp_eq*1e5/(sc.constants.c*sc.constants.h) #h/Mpc
    
    #Make the second plot.
        
    plt.figure()
    plt.plot(k, Matter_PS, label='Theory prediction')
    plt.errorbar(k_gal, P_gal, ErrorP_gal, label='SDSS Galaxies (DR7 LRG)',\
                 ls='none', marker='.', capsize=2)
    plt.errorbar(k_ACT, P_ACT, P_upper, label='CMB (WMAP+ACT)', ls='none',\
                 marker='.', capsize=2)
    plt.errorbar(k_Ly, P_Ly, DeltaP_Ly, label=r'Lyman $\alpha$ foresyt',\
                 ls='none', marker='.', capsize=2)
    plt.vlines(k_eq, plt.ylim()[0], plt.ylim()[1], label=r'$k_{eq}$',\
               ls='--')
    plt.xlabel(r'Wavenumber $k$ [$h$/Mpc]')
    plt.ylabel(r'$P(k)$ [(Mpc/$h$)$^3$]')
    plt.title('The total matter power-spectrum')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.grid()
    plt.savefig('../Plots/Milestone IV/Total matter PS.pdf')
    
def other_plots():
    
    data = np.loadtxt('../Results/Matter_PS.txt')
    
    #k: wavenumber.
    
    k: np.array(float) = data[:, 0]
    
    #keta0: k times eta_0.
    
    keta0: np.array(float) = data[:, 2]
    
    #Theta_l: photon temperature multipoles.
    
    Theta_2: np.array(float) = data[:, 3]
    
    Theta_20: np.array(float) = data[:, 4]
    
    Theta_200: np.array(float) = data[:, 5]
    
    Theta_2000: np.array(float) = data[:, 6]
    
    plt.figure()
    plt.plot(keta0, Theta_2, label=r'$\Theta_2(k)$')
    plt.plot(keta0, Theta_20, label=r'$\Theta_{20}(k)$')
    plt.plot(keta0, Theta_200, label=r'$\Theta_{200}(k)$')
    plt.plot(keta0, Theta_2000, label=r'$\Theta_{2000}(k)$')
    plt.xlabel(r'$k\eta_0$ []')
    plt.ylabel(r'$\Theta_\ell$')
    plt.title(r'$\Theta_ell$ for different $\ell$ values')
    plt.grid()
    plt.legend()
    plt.savefig('../Plots/Milestone IV/theta_l.pdf')
    
    plt.figure()
    plt.plot(k, Theta_2**2/k, label=r'$\Theta_2(k)$')
    plt.plot(k, Theta_20**2/k, label=r'$\Theta_{20}(k)$')
    plt.plot(k, Theta_200**2/k, label=r'$\Theta_{200}(k)$')
    plt.plot(k, Theta_2000**2/k, label=r'$\Theta_{2000}(k)$')
    plt.xlabel(r'$k$ []')
    plt.ylabel(r'$\Theta_\ell^2$')
    plt.title(r'$\Theta_ell^2/k$ for different $\ell$ values')
    plt.grid()
    plt.legend()
    plt.savefig('../Plots/Milestone IV/theta_l.pdf')
    
def CMB_map():
        
    nside = 128
    # Read the CMB power spectrum 'cells.txt' in folder 'Results'.
    data = np.loadtxt('../Results/cells.txt')
    # Extract the first column as l values
    ell = data[:, 0].astype(int)
    # Extract the second column as Cl values
    C_ell = data[:, 1]*(2*np.pi)/(ell*(ell+1))
    # Generate random spherical harmonic coefficients
    alm = hp.synalm(C_ell, lmax=np.max(ell), new=True)
    # Convert the spherical harmonic coefficients to a map
    cmb_map = hp.alm2map(alm, nside)
    hp.mollview(cmb_map)
    plt.show()
    
def milestone4(polarization: bool):
    
    """
    Execute Milestone IV functionality.

    This function executes Milestone IV functionality, which involves plotting
    the Cosmic Microwave Background (CMB) power-spectrum and the matter
    power-spectrum. It calls the 'CMB_PowerSpectrum' function to plot the CMB
    power-spectrum and the 'Matter_PowerSpectrum' function to plot the matter
    power-spectrum.

    Parameters:
        polarization (bool): indicates if polarization has been included.

    Returns:
        None.
    """
    
    #Plot theta_l plots.
    
    other_plots()
    
    #Plot the CMB power-spectrum.
    
    CMB_PowerSpectrum()
    
    #Plot the matter power-spectrum.
    
    Matter_PowerSpectrum(polarization)
    
    #Plot the CMB map.
    
    CMB_map()
  
milestone4(False)
