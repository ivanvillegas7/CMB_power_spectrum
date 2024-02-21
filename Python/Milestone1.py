# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 12:27:12 2024

@author: Iván Villegas Pérez
"""

import scipy as sc

import numpy as np

import matplotlib.pyplot as plt

c: float = sc.constants.c #m/s

def cosmology():
    
    data_cos = np.loadtxt('../Results/cosmology.txt')
    
    x: np.array(float) = data_cos[:, 0]
    
    eta: np.array(float) = data_cos[:, 1]
    
    Hp: np.array(float) = data_cos[:, 2]
    
    dHp: np.array(float) = data_cos[:, 3]
    
    OmegaB: np.array(float) = data_cos[:, 4]
    
    OmegaCDM: np.array(float) = data_cos[:, 5]
    
    OmegaLambda: np.array(float) = data_cos[:, 6]
    
    OmegaR: np.array(float) = data_cos[:, 7]
    
    d_L: np.array(float) = data_cos[:, 8]
    
    ddHp: np.array(float) = data_cos[:, 9]
    
    t: np.array(float) = data_cos[:, 10]
    
    z: np.array(float) = np.zeros(len(x))
    
    for i in range(len(x)):
        
        z[i] = np.exp(-x[i])-1
        
    plt.figure()
    plt.plot(x, 1/Hp*ddHp)
    plt.title(r'$\frac{1}{\mathcal{H}(x)}\cdot\frac{\text{d}^2\mathcal{H}(x}{\text{d}x^2}$ vs $x$')
    plt.grid(True)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$\frac{1}{\mathcal{H}(x)}\cdot\frac{\text{d}^2\mathcal{H}(x}{\text{d}x^2}$')
    plt.savefig('../Plots/Milestone I/ddHp_over_Hp.pdf')
    
    plt.figure()
    plt.plot(x, 1/Hp*dHp)
    plt.title(r'$\frac{1}{\mathcal{H}(x)}\cdot\frac{\text{d}\mathcal{H}(x}{\text{d}x}$ vs $x$')
    plt.grid(True)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$\frac{1}{\mathcal{H}(x)}\cdot\frac{\text{d}\mathcal{H}(x}{\text{d}x}$')
    plt.savefig('../Plots/Milestone I/dHp_over_Hp.pdf')
    
    plt.figure()
    plt.plot(x, eta*Hp/c)
    plt.title(r'$\frac{\eta(x)\mathcal{H}(x)}{c}\cdot\frac{\text{d}\mathcal{H}(x}{\text{d}x}$ vs $x$')
    plt.grid(True)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$\frac{\eta(x)\mathcal{H}(x)}{c}\cdot\frac{\text{d}\mathcal{H}(x}{\text{d}x}$')
    plt.savefig('../Plots/Milestone I/eta_times_Hp_over_c.pdf')
    
    plt.figure()
    plt.plot(x, Hp)
    plt.title(r'$\mathcal{H}(x)$ vs $x$')
    plt.grid(True)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$\mathcal{H}(x)$')
    plt.yscale('log')
    plt.savefig('../Plots/Milestone I/Hp.pdf')

    plt.figure()
    plt.plot(x, t)
    plt.title(r'$t(x)$ vs $x$')
    plt.grid(True)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$t(x)$')
    plt.savefig('../Plots/Milestone I/t.pdf')

    plt.figure()
    plt.plot(x, eta/c)
    plt.title(r'$\frac{\eta(x)}{c}$ vs $x$')
    plt.grid(True)
    plt.yscale('log')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$\frac{\eta(x)}{c}$')
    plt.savefig('../Plots/Milestone I/eta_over_c.pdf')
        
    plt.figure()
    plt.plot(x, OmegaB+OmegaCDM, label=r'$\Omega_M=\Omega_B+\Omega_{CDM}$')
    plt.plot(x, OmegaR, label=r'$\Omega_\gamma$')
    plt.plot(x, OmegaLambda, label=r'$\Omega_\Lambda$')
    plt.legend()
    plt.title(r'$\Omega_i(x)$ vs $x$')
    plt.grid(True)
    plt.xlim(-10, 0)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$\Omega_i(x)$')
    plt.savefig('../Plots/Milestone I/omegas.pdf')
    
    return(z, d_L)

def supernova():    
    
    data_sup = np.loadtxt('../data/supernovadata.txt', skiprows=1)
    
    z: np.array(float) = data_sup[:, 0]
    
    d_L: np.array(float) = data_sup[:, 1]
    
    d_L_error: np.array(float) = data_sup[:, 2]
    
    z_cos: np.array(float)
    
    d_L_cos: np.array(float)
    
    z_cos, d_L_cos = cosmology()
    
    plt.figure()
    plt.errorbar(z, d_L, d_L_error, marker='.', ls='none', label='Supernovae observations')
    plt.plot(z_cos, d_L_cos, label='Computed curve')
    plt.title(r'$d_L(z)$ vs $z$')
    plt.grid(True)
    plt.xlabel(r'$z$')
    plt.ylabel(r'$d_L(z)$')
    plt.savefig('../Plots/Milestone I/supernova_fit.pdf')
    
def MCMC_supernova_fit():    
    
    data_fit = np.loadtxt('../Results/results_supernovafitting.txt', skiprows=200)
    
    chi2: np.array(float) = data_fit[:, 0]
    
    #h: np.array(float) = data_fit[:, 1]
    
    OmegaM: np.array(float) = data_fit[:, 2]
    
    OmegaK: np.array(float) = data_fit[:, 3]
    
    chi2_min: float = min(chi2)
    
    OmegaLambda: np.array(float) = np.ones(len(OmegaM))-OmegaM-OmegaK
    
    plt.figure()
    plt.plot(OmegaM[chi2 < chi2_min + 3.53], OmegaLambda[chi2 < chi2_min + 3.53], marker='.', ls='none', label=r'$1\sigma$')
    plt.plot(np.linspace(0, 1), 1-np.linspace(0, 1), ls='dashed', label='Flat Universe', color='black')
    plt.grid(True)
    plt.xlabel(r'$\Omega_M$')
    plt.ylabel(r'$\Omega_\Lambda$')
    plt.xlim(0, 1)
    plt.ylim(0, 1.5)
    plt.savefig('../Plots/Milestone I/supernova_sigmas.pdf')
    
    plt.figure()
    plt.hist(OmegaLambda)
    plt.plot(OmegaLambda, np.random.normal(np.mean(OmegaLambda), np.std(OmegaLambda)))
    plt.xlabel(r'$\Omega_\Lambda$')
    plt.title(r'Posterior for $\Omega_\Lambda$')
    plt.vlines(OmegaLambda[chi2==chi2_min], 0, 10, linestyles='dashed', label='Best fit value')
    plt.legend()
    plt.grid(True)
    plt.savefig('../Plots/Milestone I/OmegaLambda_hist.pdf')
    
    plt.figure()
    plt.hist(OmegaM)
    plt.plot(OmegaLambda, np.random.normal(np.mean(OmegaLambda), np.std(OmegaLambda)))
    plt.xlabel(r'$\Omega_M$')
    plt.title(r'Posterior for $\Omega_M$')
    plt.vlines(OmegaM[chi2==chi2_min], 0, 10, linestyles='dashed', label='Best fit value')
    plt.legend()
    plt.grid(True)
    plt.savefig('../Plots/Milestone I/OmegaM_hist.pdf')
    
    plt.figure()
    plt.hist(OmegaK)
    plt.plot(OmegaK, np.random.normal(np.mean(OmegaK), np.std(OmegaK)))
    plt.xlabel(r'$\Omega_\K$')
    plt.title(r'Posterior for $\Omega_K$')
    plt.vlines(OmegaK[chi2==chi2_min], 0, 10, linestyles='dashed', label='Best fit value')
    plt.legend()
    plt.grid(True)
    plt.savefig('../Plots/OmegaK_hist.pdf')
       
def main():
    
    supernova()
    
    MCMC_supernova_fit()
