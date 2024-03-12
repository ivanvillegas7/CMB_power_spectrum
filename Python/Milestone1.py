# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 12:27:12 2024

@author: Iván Villegas Pérez
"""

"""
This script contains all the relevant functions for analysing and plotting the
results of running the C++ codes related to solving the background cosmology of
the Universe.

Here you can find the functions 'plot()', 'cosmology()', 'supernova()' and
'milestone1()'. All of them are explained below.
"""

#Import all relevant packages

import scipy as sc

import numpy as np

import matplotlib.pyplot as plt

#Define function 'plot()'.

def plot(x: np.array(float), y: np.array(float), i1: int, i2:int, i3: int):
    
    """
    Plot function for illustrating different epochs of cosmological evolution.

    Parameters:
        x : numpy array of float
            X-coordinates for the plot.
        y : numpy array of float
            Y-coordinates for the plot.
        i1 : int
            Index indicating the start of matter domination epoch.
        i2 : int
            Index indicating the start of dark energy domination epoch.
        i3 : int
            Index indicating the start of the accelerated expansion.

    Returns:
        None.

    This function plots different epochs of cosmological evolution based on the
    provided data.

    It fills the regions corresponding to radiation domination, matter
    domination, and dark energy domination with different colors, and marks the
    epoch limits with dashed vertical lines.
    """
    
    plt.fill_between(x[i2-1:], min(y)*np.ones(len(x)-i2+1),\
                     max(y)*np.ones(len(x)-i2+1), color='lightgreen',\
                     label=r'Dark Energy ($\Lambda$) domination')
    plt.fill_between(x[i1:i2+1], min(y)*np.ones(i2-i1+1),\
                     max(y)*np.ones(i2-i1+1), color='lightsteelblue',\
                     label='Matter domination')
    plt.fill_between(x[0:i1+1], min(y)*np.ones(i1+1), max(y)*np.ones(i1+1),\
                     color='bisque', label='Radiation domination')
    plt.vlines(x[i1], min(y), max(y), ls='dashed', color='black',\
               label='Epoch limits')
    plt.vlines(x[i2], min(y), max(y), ls='dashed', color='black')
    plt.xlim(min(x), max(x))
    plt.ylim(min(y), max(y))
    plt.legend()
            
def cosmology()->(np.array(float), np.array(float)):
    
    """
    Perform cosmological analysis and generate plots.
    
    Parameters:
        None.

    Returns:
        Tuple of numpy arrays (z, d_L).

    This function reads cosmological data from a file and performs various
    calculations and plot generations. It calculates parameters related to
    cosmological evolution such as conformal time, Hubble parameter, and
    various densities. Plots are generated for different parameters and epochs
    of cosmological evolution.

    Returns a tuple containing redshifts (z) and luminosity distances (d_L).
    """
    
    #Define the speed of light as a constant
    
    c: float = sc.constants.c #m/s
    
    #Read the data from 'cosmology.txt' in folder 'Results'.
    
    data_cos = np.loadtxt('../Results/cosmology.txt')
    
    #Assign the different variables to its corresponding data.
    
    #x=ln(a), main time variable.
    
    x: np.array(float) = data_cos[:, 0]
    
    #eta: conformal time.
    
    eta: np.array(float) = data_cos[:, 1]/(1e9*365*24*60*60)
    
    #Hp: conformal Hubble factor.
    
    Hp: np.array(float) = data_cos[:, 2]/(100*1e3/(1e6*sc.constants.parsec))
    
    #dHp: first derivative of Hp with respect to x.
    
    dHp: np.array(float) = data_cos[:, 3]/(100*1e3/(1e6*sc.constants.parsec))
    
    #OmegaB relative density of baryonic (ordinary) matter.
    
    OmegaB: np.array(float) = data_cos[:, 4]
    
    #OmegaCDM: relative density of cold dark matter.
    
    OmegaCDM: np.array(float) = data_cos[:, 5]
    
    #OmegaLambda: relative density of dark energy.
    
    OmegaLambda: np.array(float) = data_cos[:, 6]
    
    #OmegaR: relative density of radiation (photons).
    
    OmegaR: np.array(float) = data_cos[:, 7]
    
    #OmegaNu: relative density of neutrinos.

    OmegaNu: np.array(float) = data_cos[:, 8]
    
    #d_L: distante luminosity.
    
    d_L: np.array(float) = data_cos[:, 9]/(1e9*sc.constants.parsec)
    
    #ddHp: second derivative of Hp with respect to x.
    
    ddHp: np.array(float) = data_cos[:, 10]/(100*1e3/(1e6*sc.constants.parsec))
    
    #t: cosmological time.
    
    t: np.array(float) = data_cos[:, 11]/(1e9*365*24*60*60)
    
    #Define and compute redshift (z) from the given data.
    
    z: np.array(float) = np.exp(-x)-1
        
    #Define the omegas corresponding to cold and relativistic matter/particles.
        
    OmegaM:np.array(float) = OmegaB+OmegaCDM

    OmegaRel:np.array(float) = OmegaR+OmegaNu

    #Define the index for today (x=0)
    
    today: int = np.argmin(np.abs(x))
    
    #Define the second derivative of the scale factor with respect to the
    #cosmological time
    
    a_dotdot: np.array(float) = Hp**2*np.exp(-x)+\
                                (Hp[today]*np.exp(-2*x[today]))**2*\
                                (-3*OmegaM*np.exp(-x)*Hp-\
                                 4*OmegaRel*np.exp(-2*x)*Hp)/(2*Hp*np.exp(-x))
    
    #Get the indeces when the cuantity of cold matter is the same (or most 
    #similar) to the cuantity of relativistic particles.
    
    index_M_R: int = np.argmin(np.abs(OmegaRel-OmegaM))
    
    index_M_Lambda: int = np.argmin(np.abs(OmegaLambda-OmegaM))
    
    #Get the index when the Universe starts its accelerated expansion.
    
    index: int = np.argmin(np.abs(a_dotdot))
    
    #Make the different plots and save them. 
        
    plt.figure()
    plt.plot(x, ddHp/Hp)
    plot(x, ddHp/Hp, index_M_R, index_M_Lambda, index)
    plt.title(r'$\frac{1}{\mathcal{H}(x)}\cdot\frac{\text{d}^2\mathcal{H}(x)}{\text{d}x^2}$ vs $x$')
    plt.grid(True)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$\frac{1}{\mathcal{H}(x)}\cdot\frac{\text{d}^2\mathcal{H}(x)}{\text{d}x^2}$')
    plt.savefig('../Plots/Milestone I/ddHp_over_Hp.pdf')
    
    plt.figure()
    plt.plot(x, dHp/Hp)
    plot(x, dHp/Hp, index_M_R, index_M_Lambda, index)
    plt.title(r'$\frac{1}{\mathcal{H}(x)}\cdot\frac{\text{d}\mathcal{H}(x)}{\text{d}x}$ vs $x$')
    plt.grid(True)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$\frac{1}{\mathcal{H}(x)}\cdot\frac{\text{d}\mathcal{H}(x)}{\text{d}x}$')
    plt.savefig('../Plots/Milestone I/dHp_over_Hp.pdf')
        
    plt.figure()
    plt.plot(x, eta*Hp/(c*1.02*1e-1))
    plot(x, eta*Hp/(c*1.02*1e-1), index_M_R, index_M_Lambda, index)
    plt.title(r'$100\frac{\eta(x)\mathcal{H}(x)}{c}$ vs $x$')
    plt.grid(True)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$100\frac{\eta(x)\mathcal{H}(x)}{c}$')
    plt.savefig('../Plots/Milestone I/eta_times_Hp_over_c.pdf')
    
    plt.figure()
    plt.plot(x, 100*Hp)
    plot(x, 100*Hp, index_M_R, index_M_Lambda, index)
    plt.title(r'$\mathcal{H}(x)$ vs $x$')
    plt.grid(True)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$\mathcal{H}(x)$ [$\frac{\text{ km}}{\text{Mpc s}}$]')
    plt.yscale('log')
    plt.savefig('../Plots/Milestone I/Hp.pdf')

    plt.figure()
    plt.plot(x, t)
    plot(x, t, index_M_R, index_M_Lambda, index)
    plt.title(r'$t(x)$ vs $x$')
    plt.grid(True)
    plt.yscale('log')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$t(x)$ [Gyr]')
    plt.savefig('../Plots/Milestone I/t.pdf')

    plt.figure()
    plt.plot(x, eta/c)
    plot(x, eta/c, index_M_R, index_M_Lambda, index)
    plt.title(r'$\frac{\eta(x)}{c}$ vs $x$')
    plt.grid(True)
    plt.yscale('log')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$\frac{\eta(x)}{c}$ [Gyr]')
    plt.savefig('../Plots/Milestone I/eta_over_c.pdf')
    
    plt.figure()    
    plt.plot(x, OmegaLambda, label=r'$\Omega_\Lambda$', color='green')
    plt.fill_between(x[index_M_Lambda-1:], np.zeros(len(x)-index_M_Lambda+1),\
                     np.ones(len(x)-index_M_Lambda+1), color='lightgreen',\
                     label=r'Dark Energy ($\Lambda$) domination')
    plt.plot(x, OmegaM, label=r'$\Omega_M=\Omega_B+\Omega_{CDM}$',\
             color='royalblue')
    plt.plot(x, OmegaB, label=r'$\Omega_B$', ls='dashed', color='darkblue')
    plt.plot(x, OmegaCDM, label=r'$\Omega_{CDM}$', ls='dashed', color='cyan')
    plt.fill_between(x[index_M_R:index_M_Lambda+1],\
                     np.zeros(index_M_Lambda-index_M_R+1),\
                     np.ones(index_M_Lambda-index_M_R+1),\
                     color='lightsteelblue', label='Matter domination')
    plt.plot(x, OmegaRel, label=r'$\Omega_R=\Omega_\gamma+\Omega_\nu$',\
             color='darkorange')
    plt.plot(x, OmegaNu, label=r'$\Omega_\nu$', ls='dashed', color='red')
    plt.plot(x, OmegaR, label=r'$\Omega_\gamma$', ls='dashed', color='y')
    plt.fill_between(x[0:index_M_R+1], np.zeros(index_M_R+1),\
                     np.ones(index_M_R+1), color='bisque',\
                     label='Radiation domination')
    plt.legend()
    plt.title(r'$\Omega_i(x)$ vs $x$')
    plt.grid(True)
    plt.xlim(-15, 0)
    plt.ylim(0, 1)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$\Omega_i(x)$')
    plt.savefig('../Plots/Milestone I/omegas.pdf')
    
    #Print the times (x, redshift and cosmological time) when there is a change
    #of domination. Print the age of the Universe and the conformal time today.
    
    print(f'\nTIMES\n\
                  x       z       t [Gyr]\n\
          M-R:  {x[index_M_R]:.2f}  {z[index_M_R]:.2f}     {t[index_M_R]*1e5:.2f}e-5\n\
          M_Λ:  {x[index_M_Lambda]:.2f}     {z[index_M_Lambda]:.2f}    {t[index_M_Lambda]:.2f}\n\
          ä=0:  {x[index]:.2f}     {z[index]:.2f}     {t[index]:.2f}\n')
    
    print(f'Age of the Universe: t(0)≈{t[today]:.2f} Gyr.\n')
    
    print(f'Conformal time today: η(0)/c≈{eta[today]/c:.2f} Gyr.\n')
    
    #Return redshift (z) and luminosity distance (d_L) in a tuple.
    
    return(z, d_L)

def MCMC_supernova_fit():
    
    """
    Perform analysis and generate plots for supernova fitting results obtained
    via MCMC.

    This function reads results from a file, performs analysis on parameters
    OmegaM, OmegaLambda, and OmegaK, and generates histograms and scatter plots
    to visualize the posterior distributions.
    
    Parameters:
        None.

    Returns:
        None.
    """
    
    #Read the data from 'results_supernovafitting.txt' in folder 'Results'.
    # Ignore first 200 rows.
    
    data = np.loadtxt('../Results/results_supernovafitting.txt', skiprows=200)
    
    #Assign the different variables to its corresponding data.
    
    #chi2: χ² values for each combination of fitting parameters
    
    chi2: np.array(float) = data[:, 0]
    
    #H: fitted Hubble factor parameter.
    
    H: np.array(float) = 100*data[:, 1]
    
    #OmegaM: fitted cold matter parameter.
    
    OmegaM: np.array(float) = data[:, 2]
    
    #OmegaK: fitted curvature parameter.
    
    OmegaK: np.array(float) = data[:, 3]
    
    #chi2_min: minimum value of χ², identifying the index of the best fit.
    
    chi2_min: float = min(chi2)
    
    index: int = np.argmin(chi2)
    
    #Compute the corresponding fitting parameter for dark energy (OmegaLambda).
    
    OmegaLambda: np.array(float) = np.ones(len(OmegaM))-OmegaM-OmegaK
    
    #Plot the 1σ and 2σ deviation from the best fit, along with rhe values 
    #compatible with a flat Universe.
    
    plt.figure()
    plt.plot(OmegaM[chi2 < chi2_min + 8.02],\
             OmegaLambda[chi2 < chi2_min + 8.02], marker='.', ls='none',\
                 label=r'$2\sigma$')
    plt.plot(OmegaM[chi2 < chi2_min + 3.53],\
             OmegaLambda[chi2 < chi2_min + 3.53], marker='.', ls='none',\
                 label=r'$1\sigma$')
    plt.plot(np.linspace(0, 1), 1-np.linspace(0, 1), ls='dashed',\
             label='Flat Universe', color='black')
    plt.grid(True)
    plt.legend()
    plt.xlabel(r'$\Omega_M$')
    plt.ylabel(r'$\Omega_\Lambda$')
    plt.xlim(0, 1)
    plt.ylim(0, 1.4)
    plt.savefig('../Plots/Milestone I/supernova_sigmas.pdf')
    
    #Set the number of bins for the histograms.
    
    n_bins: int = 15
    
    #Define the parameters for the gaussian fit for OmegaLambda, compute them
    #and get the gaussian distribution.
    
    mu_Lambda: float
    
    std_Lambda: float
    
    mu_Lambda, std_Lambda = sc.stats.norm.fit(OmegaLambda)
    
    pdf_Lambda: np.array(float) = sc.stats.norm.pdf(OmegaLambda, mu_Lambda,\
                                                    std_Lambda)
        
    #Plot and save the histogram.
    
    plt.figure()
    plt.hist(OmegaLambda, bins=n_bins, density=True)
    plt.plot(OmegaLambda, pdf_Lambda, marker='o', ls='none')
    plt.xlabel(r'$\Omega_\Lambda$')
    plt.title(r'Posterior for $\Omega_\Lambda$')
    plt.vlines(OmegaLambda[index], 0, max(pdf_Lambda),\
               linestyles='dashed', label='Best fit value', color='black')
    plt.legend()
    plt.grid(True)
    plt.savefig('../Plots/Milestone I/OmegaLambda_hist.pdf')
    
    #Define the parameters for the gaussian fit for OmegaM, compute them
    #and get the gaussian distribution.
    
    mu_M: float
    
    std_M: float
    
    mu_M, std_M = sc.stats.norm.fit(OmegaM)
    
    pdf_M: np.array(float) = sc.stats.norm.pdf(OmegaM, mu_M, std_M)
    
    #Plot and save the histogram.
    
    plt.figure()
    plt.hist(OmegaM, bins=n_bins, density=True)
    plt.plot(OmegaM, pdf_M, marker='o', ls='none')
    plt.xlabel(r'$\Omega_M$')
    plt.title(r'Posterior for $\Omega_M$')
    plt.vlines(OmegaM[index], 0, max(pdf_M), linestyles='dashed',\
               label='Best fit value', color='black')
    plt.legend()
    plt.grid(True)
    plt.savefig('../Plots/Milestone I/OmegaM_hist.pdf')
    
    #Define the parameters for the gaussian fit for OmegaK, compute them
    #and get the gaussian distribution.
    
    mu_K: float
    
    std_K: float
    
    mu_K, std_K = sc.stats.norm.fit(OmegaK)
    
    pdf_K: np.array(float) = sc.stats.norm.pdf(OmegaK, mu_K, std_K)
    
    #Plot and save the histogram.
    
    plt.figure()
    plt.hist(OmegaK, bins=n_bins, density=True)
    plt.plot(OmegaK, pdf_K, marker='o', ls='none')
    plt.xlabel(r'$\Omega_K$')
    plt.title(r'Posterior for $\Omega_K$')
    plt.vlines(OmegaK[index], 0, max(pdf_K), linestyles='dashed',\
               label='Best fit value', color='black')
    plt.legend()
    plt.grid(True)
    plt.savefig('../Plots/Milestone I/OmegaK_hist.pdf')
    
    #Define the parameters for the gaussian fit for OmegaK, compute them
    #and get the gaussian distribution.
    
    mu_H: float
    
    std_H: float
    
    mu_H, std_H = sc.stats.norm.fit(H)
    
    pdf_H: np.array(float) = sc.stats.norm.pdf(H, mu_H, std_H)
    
    #Plot and save the histogram.
    
    plt.figure()
    plt.hist(H, bins=n_bins, density=True)
    plt.plot(H, pdf_H, marker='o', ls='none')
    plt.xlabel(r'$H_0$')
    plt.title(r'Posterior for $H_0$')
    plt.vlines(H[index], 0, max(pdf_H), linestyles='dashed',\
               label='Best fit value', color='black')
    plt.legend()
    plt.grid(True)
    plt.savefig('../Plots/Milestone I/H0_hist.pdf')

    #Print thge best fitting parameters

    print(f'\nPARAMETERS\n\
          χ²     Ω_M     Ω_k     Ω_Λ      H_0\n\
        {chi2_min:.2f}    {OmegaM[index]:.2f}    {OmegaK[index]:.2f}    {OmegaLambda[index]:.2f}    {H[index]:.2f}\n')
    
def supernova():
    
    """
    Perform analysis and generate plots for supernova data fitting.

    This function reads supernova data from a file, computes the luminosity
    distance using cosmological calculations, and generates a plot comparing
    the computed curve with the observed supernova data.
    
    Parameters:
        None.

    Returns:
        None.
    """
    
    #Read the data from 'supernovadata.txt' in folder 'data'. Ignore first row.
    
    data_sup = np.loadtxt('../data/supernovadata.txt', skiprows=1)
    
    #Assign the different variables to its corresponding data.
    
    #z: redshift
    
    z: np.array(float) = data_sup[:, 0]
    
    #d_L: luminosity distance.
    
    d_L: np.array(float) = data_sup[:, 1]
    
    #d_L_error: luminosity distance error.
    
    d_L_error: np.array(float) = data_sup[:, 2]
    
    #Define and get both redshift and lumnosity distance from the function
    #'cosmology()'.
    
    z_cos: np.array(float)
    
    d_L_cos: np.array(float)
    
    z_cos, d_L_cos = cosmology()
    
    #Plot the observational data (with error bars) and the fit from the data
    #given by the function 'cosmology()'.
    
    plt.figure()
    plt.errorbar(z, d_L, d_L_error, marker='.', ls='none',\
                 label='Supernovae observations')
    plt.plot(z_cos, d_L_cos, label='Fiduicial cosmology curve')
    plt.title(r'$d_L(z)$ vs $z$')
    plt.legend()
    plt.grid(True)
    plt.xlabel(r'$z$')
    plt.ylabel(r'$d_L(z)$ [Gpc]')
    plt.xlim(0, 1.4)
    plt.ylim(0, 10)
    plt.savefig('../Plots/Milestone I/supernova_fit.pdf')
    
    plt.figure()
    plt.errorbar(z, d_L/z, d_L_error/z, marker='.', ls='none',\
                 label='Supernovae observations')
    plt.plot(z_cos[z_cos>0.015], d_L_cos[z_cos>0.015]/z_cos[z_cos>0.015],\
             label='Fiduicial cosmology curve')
    plt.title(r'$d_L(z)/z$ vs $z$')
    plt.legend()
    plt.grid(True)
    plt.xlabel(r'$z$')
    plt.ylabel(r'$d_L(z)/z$ [Gpc]')
    plt.xlim(0, 1.4)
    plt.ylim(3.5, 8)
    plt.savefig('../Plots/Milestone I/supernova_fit_over_z.pdf')
       
def milestone1():
    
    """
    Execute main functionality for supernova analysis and background cosmology
    solving.

    This function serves as the entry point for the supernova analysis program.
    It calls the 'supernova()' function to analyze supernova data and generate
    plots, followed by the 'MCMC_supernova_fit()' function to show the best
    fitting parameters computed after having used Markov Chain Monte Carlo.
    
    Parameters:
        None.

    Returns:
        None.
    """
    
    #Run the functions.
    
    MCMC_supernova_fit()
    
    supernova()
