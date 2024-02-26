# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 12:27:12 2024

@author: Iván Villegas Pérez
"""

import scipy as sc

import numpy as np

import matplotlib.pyplot as plt

def plot(x: np.array(float), y: np.array(float), i1: int, i2:int):
    
    plt.fill_between(x[i2-1:], min(y)*np.ones(len(x)-i2+1),\
                     max(y)*np.ones(len(x)-i2+1), color='lightgreen',\
                     label=r'Dark Eenrgy ($\Lambda$) domination')
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
    
    c: float = sc.constants.c #m/s
    
    data_cos = np.loadtxt('../Results/cosmology.txt')
    
    x: np.array(float) = data_cos[:, 0]
    
    eta: np.array(float) = data_cos[:, 1]/(1e9*365*24*60*60)
    
    Hp: np.array(float) = data_cos[:, 2]/(100*1e3/(1e9*sc.constants.parsec))
    
    dHp: np.array(float) = data_cos[:, 3]/(100*1e3/(1e9*sc.constants.parsec))
    
    OmegaB: np.array(float) = data_cos[:, 4]
    
    OmegaCDM: np.array(float) = data_cos[:, 5]
    
    OmegaLambda: np.array(float) = data_cos[:, 6]
    
    OmegaR: np.array(float) = data_cos[:, 7]

    OmegaNu: np.array(float) = data_cos[:, 8]
    
    d_L: np.array(float) = data_cos[:, 9]/(1e9*sc.constants.parsec)
    
    ddHp: np.array(float) = data_cos[:, 10]/(100*1e3/(1e9*sc.constants.parsec))
    
    t: np.array(float) = data_cos[:, 11]/(1e9*365*24*60*60)
    
    z: np.array(float) = np.zeros(len(x))
    
    for i in range(len(x)):
        
        z[i] = np.exp(-x[i])-1
        
    OmegaM:np.array(float) = OmegaB+OmegaCDM

    OmegaRel:np.array(float) = OmegaR+OmegaNu
    
    index_M_R: int = np.argmin(np.abs(OmegaRel-OmegaM))
    
    index_M_Lambda: int = np.argmin(np.abs(OmegaLambda-OmegaM))
        
    plt.figure()
    plt.plot(x, ddHp/Hp)
    plot(x, ddHp/Hp, index_M_R, index_M_Lambda)
    plt.title(r'$\frac{1}{\mathcal{H}(x)}\cdot\frac{\text{d}^2\mathcal{H}(x)}{\text{d}x^2}$ vs $x$')
    plt.grid(True)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$\frac{1}{\mathcal{H}(x)}\cdot\frac{\text{d}^2\mathcal{H}(x)}{\text{d}x^2}$')
    plt.savefig('../Plots/Milestone I/ddHp_over_Hp.pdf')
    
    plt.figure()
    plt.plot(x, dHp/Hp)
    plot(x, dHp/Hp, index_M_R, index_M_Lambda)
    plt.title(r'$\frac{1}{\mathcal{H}(x)}\cdot\frac{\text{d}\mathcal{H}(x)}{\text{d}x}$ vs $x$')
    plt.grid(True)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$\frac{1}{\mathcal{H}(x)}\cdot\frac{\text{d}\mathcal{H}(x)}{\text{d}x}$')
    plt.savefig('../Plots/Milestone I/dHp_over_Hp.pdf')
    
    plt.figure()
    plt.plot(x, eta*Hp/c)
    plot(x, eta*Hp/c, index_M_R, index_M_Lambda)
    plt.title(r'$\frac{\eta(x)\mathcal{H}(x)}{c}\cdot\frac{\text{d}\mathcal{H}(x)}{\text{d}x}$ vs $x$')
    plt.grid(True)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$\frac{\eta(x)\mathcal{H}(x)}{c}\cdot\frac{\text{d}\mathcal{H}(x)}{\text{d}x}$')
    plt.savefig('../Plots/Milestone I/eta_times_Hp_over_c.pdf')
    
    plt.figure()
    plt.plot(x, Hp)
    plot(x, Hp, index_M_R, index_M_Lambda)
    plt.title(r'$\mathcal{H}(x)$ vs $x$')
    plt.grid(True)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$\mathcal{H}(x)$')
    plt.yscale('log')
    plt.savefig('../Plots/Milestone I/Hp.pdf')

    plt.figure()
    plt.plot(x, t)
    plot(x, t, index_M_R, index_M_Lambda)
    plt.title(r'$t(x)$ vs $x$')
    plt.grid(True)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$t(x)$')
    plt.savefig('../Plots/Milestone I/t.pdf')

    plt.figure()
    plt.plot(x, eta/c)
    plot(x, eta/c, index_M_R, index_M_Lambda)
    plt.title(r'$\frac{\eta(x)}{c}$ vs $x$')
    plt.grid(True)
    plt.yscale('log')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$\frac{\eta(x)}{c}$')
    plt.savefig('../Plots/Milestone I/eta_over_c.pdf')
    
    plt.figure()    
    plt.plot(x, OmegaLambda, label=r'$\Omega_\Lambda$', color='green')
    plt.fill_between(x[index_M_Lambda-1:], np.zeros(len(x)-index_M_Lambda+1),\
                     1*np.ones(len(x)-index_M_Lambda+1), color='lightgreen',\
                     label=r'Dark Energy ($\Lambda$) domination')
    plt.plot(x, OmegaM, label=r'$\Omega_M=\Omega_B+\Omega_{CDM}$',\
             color='royalblue')
    plt.fill_between(x[index_M_R:index_M_Lambda+1],\
                     np.zeros(index_M_Lambda-index_M_R+1),\
                     1*np.ones(index_M_Lambda-index_M_R+1),\
                     color='lightsteelblue', label='Matter domination')
    plt.plot(x, OmegaRel, label=r'$\Omega_R=\Omega_\gamma+\Omega_\nu$',\
             color='darkorange')
    plt.fill_between(x[0:index_M_R+1], np.zeros(index_M_R+1),\
                     1*np.ones(index_M_R+1), color='bisque',\
                     label='Radiation domination')
    plt.legend()
    plt.title(r'$\Omega_i(x)$ vs $x$')
    plt.grid(True)
    plt.xlim(-10, 0)
    plt.ylim(0, 1)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$\Omega_i(x)$')
    plt.savefig('../Plots/Milestone I/omegas.pdf')
    
    print(f'\nTIMES\n\
                  x     z         t [Gyr]\n\
          M-R:  {x[index_M_R]:.2f}  {z[index_M_R]:.2f}     {t[index_M_R]:.2f}\n\
          M_Λ:   {x[index_M_Lambda]:.2f}    {z[index_M_Lambda]:.2f}    {t[index_M_Lambda]:.2f}\n')
    
    print(f'Age of the Universe: t(0)≈{t[-1]:.2f} Gyr.\n')
    
    print(f'Conformal time today: η(0)/c≈{eta[-1]/c:.2f} Gyr.\n')
    
    return(z, d_L)

def MCMC_supernova_fit():
    
    data = np.loadtxt('../Results/results_supernovafitting.txt', skiprows=200)
    
    chi2: np.array(float) = data[:, 0]
    
    OmegaM: np.array(float) = data[:, 2]
    
    OmegaK: np.array(float) = data[:, 3]
    
    chi2_min: float = min(chi2)
    
    OmegaLambda: np.array(float) = np.ones(len(OmegaM))-OmegaM-OmegaK
    
    plt.figure()
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
    
    n_bins: int = 15
    
    mu_Lambda: float
    
    std_Lambda: float
    
    mu_Lambda, std_Lambda = sc.stats.norm.fit(OmegaLambda)
    
    pdf_Lambda: np.array(float) = sc.stats.norm.pdf(OmegaLambda, mu_Lambda,\
                                                    std_Lambda)
    
    plt.figure()
    plt.hist(OmegaLambda, bins=n_bins, density=True)
    plt.plot(OmegaLambda, pdf_Lambda, marker='o', ls='none')
    plt.xlabel(r'$\Omega_\Lambda$')
    plt.title(r'Posterior for $\Omega_\Lambda$')
    plt.vlines(OmegaLambda[chi2==chi2_min], 0, max(pdf_Lambda),\
               linestyles='dashed', label='Best fit value', color='black')
    plt.legend()
    plt.grid(True)
    plt.savefig('../Plots/Milestone I/OmegaLambda_hist.pdf')
    
    mu_M: float
    
    std_M: float
    
    mu_M, std_M = sc.stats.norm.fit(OmegaM)
    
    pdf_M: np.array(float) = sc.stats.norm.pdf(OmegaM, mu_M, std_M)
    
    plt.figure()
    plt.hist(OmegaM, bins=n_bins, density=True)
    plt.plot(OmegaM, pdf_M, marker='o', ls='none')
    plt.xlabel(r'$\Omega_M$')
    plt.title(r'Posterior for $\Omega_M$')
    plt.vlines(OmegaM[chi2==chi2_min], 0, max(pdf_M), linestyles='dashed',\
               label='Best fit value', color='black')
    plt.legend()
    plt.grid(True)
    plt.savefig('../Plots/Milestone I/OmegaM_hist.pdf')
    
    mu_K: float
    
    std_K: float
    
    mu_K, std_K = sc.stats.norm.fit(OmegaK)
    
    pdf_K: np.array(float) = sc.stats.norm.pdf(OmegaK, mu_K, std_K)
    
    plt.figure()
    plt.hist(OmegaK, bins=n_bins, density=True)
    plt.plot(OmegaK, pdf_K, marker='o', ls='none')
    plt.xlabel(r'$\Omega_K$')
    plt.title(r'Posterior for $\Omega_K$')
    plt.vlines(OmegaK[chi2==chi2_min], 0, max(pdf_K), linestyles='dashed',\
               label='Best fit value', color='black')
    plt.legend()
    plt.grid(True)
    plt.savefig('../Plots/Milestone I/OmegaK_hist.pdf')
    
def supernova():    
    
    data_sup = np.loadtxt('../data/supernovadata.txt', skiprows=1)
    
    z: np.array(float) = data_sup[:, 0]
    
    d_L: np.array(float) = data_sup[:, 1]
    
    d_L_error: np.array(float) = data_sup[:, 2]
    
    z_cos: np.array(float)
    
    d_L_cos: np.array(float)
    
    z_cos, d_L_cos = cosmology()
    
    plt.figure()
    plt.errorbar(z, d_L, d_L_error, marker='.', ls='none',\
                 label='Supernovae observations')
    plt.plot(z_cos, d_L_cos, label='Computed curve')
    plt.title(r'$d_L(z)$ vs $z$')
    plt.legend()
    plt.grid(True)
    plt.xlabel(r'$z$')
    plt.ylabel(r'$d_L(z)$')
    plt.xlim(0, 1.4)
    plt.ylim(0, 10)
    plt.savefig('../Plots/Milestone I/supernova_fit.pdf')
       
def main():
    
    supernova()
    
    MCMC_supernova_fit()
