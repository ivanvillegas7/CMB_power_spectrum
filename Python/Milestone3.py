# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 12:06:38 2024

@author: Iván Villegas Pérez
"""

#Import all relevant packages

import numpy as np

import matplotlib.pyplot as plt

def plot():
    
    polarization_txt: str
    
    polarization_txt = input('Did you include photon polarization? [Yes, no]: ')
    
    polarization: bool
    
    if polarization_txt.lower()=='no':
        
        polarization = False
        
    else:
        
        polarization = True
        
    neutrinos_txt: str
    
    neutrinos_txt = input('Did you include neutrinos? [Yes, no]: ')
    
    neutrinos: bool
    
    if neutrinos_txt.lower()=='no':
        
        neutrinos = False
        
    else:
        
        neutrinos = True
    
    #Read the data from 'perturbations_k1.txt' in folder 'Results'.
    
    data_1 = np.loadtxt('../Results/perturbations_k1.txt')
    
    #x=ln(a): main time variable.
    
    x: np.array(float) = data_1[:, 0]
    
    Theta_0_1: np.array(float) = data_1[:, 1]
    Theta_1_1: np.array(float) = data_1[:, 2]
    Theta_2_1: np.array(float) = data_1[:, 3]
    Phi_1: np.array(float) = data_1[:, 4]
    Psi_1: np.array(float) = data_1[:, 5]
    Pi_1: np.array(float) = data_1[:, 6]
    Source_T_0_1: np.array(float) = data_1[:, 7]
    Source_T_1_1: np.array(float) = data_1[:, 8]
    Source_T_2_1: np.array(float) = data_1[:, 9]
    Source_T_3_1: np.array(float) = data_1[:, 10]
    
    #Read the data from 'perturbations_k01.txt' in folder 'Results'.
    
    data_01 = np.loadtxt('../Results/perturbations_k01.txt')
    
    Theta_0_01: np.array(float) = data_01[:, 1]
    Theta_1_01: np.array(float) = data_01[:, 2]
    Theta_2_01: np.array(float) = data_01[:, 3]
    Phi_01: np.array(float) = data_01[:, 4]
    Psi_01: np.array(float) = data_01[:, 5]
    Pi_01: np.array(float) = data_01[:, 6]
    Source_T_0_01: np.array(float) = data_01[:, 7]
    Source_T_1_01: np.array(float) = data_01[:, 8]
    Source_T_2_01: np.array(float) = data_01[:, 9]
    Source_T_3_01: np.array(float) = data_01[:, 10]
    
    #Read the data from 'perturbations_k001.txt' in folder 'Results'.
    
    data_001 = np.loadtxt('../Results/perturbations_k001.txt')
    
    Theta_0_001: np.array(float) = data_001[:, 1]
    Theta_1_001: np.array(float) = data_001[:, 2]
    Theta_2_001: np.array(float) = data_001[:, 3]
    Phi_001: np.array(float) = data_001[:, 4]
    Psi_001: np.array(float) = data_001[:, 5]
    Pi_001: np.array(float) = data_001[:, 6]
    Source_T_0_001: np.array(float) = data_001[:, 7]
    Source_T_1_001: np.array(float) = data_001[:, 8]
    Source_T_2_001: np.array(float) = data_001[:, 9]
    Source_T_3_001: np.array(float) = data_001[:, 10]
    
    #Read the data from 'perturbations_k0001.txt' in folder 'Results'.
    
    data_0001 = np.loadtxt('../Results/perturbations_k0001.txt')
    
    Theta_0_0001: np.array(float) = data_0001[:, 1]
    Theta_1_0001: np.array(float) = data_0001[:, 2]
    Theta_2_0001: np.array(float) = data_0001[:, 3]
    Phi_0001: np.array(float) = data_0001[:, 4]
    Psi_0001: np.array(float) = data_0001[:, 5]
    Pi_0001: np.array(float) = data_0001[:, 6]
    Source_T_0_0001: np.array(float) = data_0001[:, 7]
    Source_T_1_0001: np.array(float) = data_0001[:, 8]
    Source_T_2_0001: np.array(float) = data_0001[:, 9]
    Source_T_3_0001: np.array(float) = data_0001[:, 10]
    
    plt.figure()
    plt.plot(x, 4*Theta_0_1, label=r'$\delta_\gamma(k=0.1/\text{Mpc})$')
    plt.plot(x, 4*Theta_0_01, label=r'$\delta_\gamma(k=0.01/\text{Mpc})$')
    plt.plot(x, 4*Theta_0_001, label=r'$\delta_\gamma(k=0.001/\text{Mpc})$')
    plt.plot(x, 4*Theta_0_0001, label=r'$\delta_\gamma(k=0.0001/\text{Mpc})$')
    plt.plot(x, delta_CDM_1, label=r'$\delta_\text{CDM}(k=0.1/\text{Mpc})$')
    plt.plot(x, delta_CDM_01, label=r'$\delta_\text{CDM}(k=0.01/\text{Mpc})$')
    plt.plot(x, delta_CDM_001, label=r'$\delta_\text{CDM}(k=0.001/\text{Mpc})$')
    plt.plot(x, delta_CDM_0001, label=r'$\delta_\text{CDM}(k=0.0001/\text{Mpc})$')
    plt.plot(x, delta_B_1, label=r'$\delta_\text{B}(k=0.1/\text{Mpc})$')
    plt.plot(x, delta_B_01, label=r'$\delta_\text{B}(k=0.01/\text{Mpc})$')
    plt.plot(x, delta_B_001, label=r'$\delta_\text{B}(k=0.001/\text{Mpc})$')
    plt.plot(x, delta_B_0001, label=r'$\delta_\text{B}(k=0.0001/\text{Mpc})$')
    if neutrinos:
        plt.plot(x, 4*Nu_0_1, label=r'$\delta_\nu(k=0.1/\text{Mpc})$')
        plt.plot(x, 4*Nu_0_01, label=r'$\delta_\nu(k=0.01/\text{Mpc})$')
        plt.plot(x, 4*Nu_0_001, label=r'$\delta_\nu(k=0.001/\text{Mpc})$')
        plt.plot(x, 4*Nu_0_0001, label=r'$\delta_\nu(k=0.0001/\text{Mpc})$')
    plt.xlabel(r'$x$')
    plt.ylabel(r'Density perturbation ($\delta_i$)')
    plt.title(r'Evolution of density perturbations over time ($x$)')
    plt.legend()
    plt.grid(True)
    plt.savefig('../Plots/Milestone III/density perturbations.pdf')
    
    plt.figure()
    plt.plot(x, -3*Theta_1_1, label=r'$v_\gamma(k=0.1/\text{Mpc})$')
    plt.plot(x, -3*Theta_1_01, label=r'$v_\gamma(k=0.01/\text{Mpc})$')
    plt.plot(x, -3*Theta_1_001, label=r'$v_\gamma(k=0.001/\text{Mpc})$')
    plt.plot(x, -3*Theta_1_0001, label=r'$v_\gamma(k=0.0001/\text{Mpc})$')
    plt.plot(x, v_CDM_1, label=r'$v_\text{CDM}(k=0.1/\text{Mpc})$')
    plt.plot(x, v_CDM_01, label=r'$v_\text{CDM}(k=0.01/\text{Mpc})$')
    plt.plot(x, v_CDM_001, label=r'$v_\text{CDM}(k=0.001/\text{Mpc})$')
    plt.plot(x, v_CDM_0001, label=r'$v_\text{CDM}(k=0.0001/\text{Mpc})$')
    plt.plot(x, v_B_1, label=r'$v_\text{B}(k=0.1/\text{Mpc})$')
    plt.plot(x, v_B_01, label=r'$v_\text{B}(k=0.01/\text{Mpc})$')
    plt.plot(x, v_B_001, label=r'$v_\text{B}(k=0.001/\text{Mpc})$')
    plt.plot(x, v_B_0001, label=r'$v_\text{B}(k=0.0001/\text{Mpc})$')
    if neutrinos:
        plt.plot(x, -3*Nu_1_1, label=r'$v_\nu(k=0.1/\text{Mpc})$')
        plt.plot(x, -3*Nu_1_01, label=r'$v_\nu(k=0.01/\text{Mpc})$')
        plt.plot(x, -3*Nu_1_001, label=r'$v_\nu(k=0.001/\text{Mpc})$')
        plt.plot(x, -3*Nu_1_0001, label=r'$v_\nu(k=0.0001/\text{Mpc})$')
    plt.xlabel(r'$x$')
    plt.ylabel(r'Velocity perturbation ($v_i$)')
    plt.title(r'Evolution of velocity perturbations over time ($x$)')
    plt.legend()
    plt.grid(True)
    plt.savefig('../Plots/Milestone III/velocity perturbations.pdf')
    
    plt.figure()
    plt.plot(x, Theta_2_1, label=r'$\Theta_2(k=0.1/\text{Mpc})$')
    plt.plot(x, Theta_2_01, label=r'$\Theta_2(k=0.01/\text{Mpc})$')
    plt.plot(x, Theta_2_001, label=r'$\Theta_2(k=0.001/\text{Mpc})$')
    plt.plot(x, Theta_2_0001, label=r'$\Theta_2(k=0.0001/\text{Mpc})$')
    if neutrinos:
        plt.plot(x, Nu_2_1, label=r'$\mathcal{N}_2(k=0.1/\text{Mpc})$')
        plt.plot(x, Nu_2_01, label=r'$\mathcal{N}_2(k=0.01/\text{Mpc})$')
        plt.plot(x, Nu_2_001, label=r'$\mathcal{N}_2(k=0.001/\text{Mpc})$')
        plt.plot(x, Nu_2_0001, label=r'$\mathcal{N}_2(k=0.0001/\text{Mpc})$')
    plt.xlabel(r'$x$')
    plt.ylabel(r'Quadrupoles')
    plt.title(r'Evolution of quadrupoles over time ($x$)')
    plt.legend()
    plt.grid(True)
    plt.savefig('../Plots/Milestone III/quadrupoles.pdf')
    
    plt.figure()
    plt.plot(x, Phi_1, label=r'$\Phi(k=0.1/\text{Mpc})$')
    plt.plot(x, Phi_1+Psi_1,\
             label=r'$\Phi(k=0.1/\text{Mpc})+\Psi(k=0.1/\text{Mpc})$')
    plt.plot(x, Phi_01, label=r'$\Phi(k=0.01/\text{Mpc})$')
    plt.plot(x, Phi_01+Psi_01,\
             label=r'$\Phi(k=0.01/\text{Mpc})+\Psi(k=0.01/\text{Mpc})$')
    plt.plot(x, Phi_001,\
             label=r'$\Phi(k=0.001/\text{Mpc})$')
    plt.plot(x, Phi_001+Psi_001,\
             label=r'$\Phi(k=0.001/\text{Mpc})+\Psi(k=0.001/\text{Mpc})$')
    plt.plot(x, Phi_0001,\
             label=r'$\Phi(k=0.0001/\text{Mpc})$')
    plt.plot(x, Phi_0001+Psi_0001,\
             label=r'$\Phi(k=0.0001/\text{Mpc})+\Psi(k=0.0001/\text{Mpc})$')
    plt.xlabel(r'$x$')
    plt.ylabel(r'Potentials')
    plt.title(r'Evolution of potentials over time ($x$)')
    plt.legend()
    plt.grid(True)
    plt.savefig('../Plots/Milestone III/potentials.pdf')
        
    if polarization:
        
        plt.figure()
        plt.plot(x, ThetaP_0_1, label=r'$\Theta_0^P(k=0.1/\text{Mpc})$')
        plt.plot(x, ThetaP_0_01, label=r'$\Theta_0^P(k=0.01/\text{Mpc})$')
        plt.plot(x, ThetaP_0_001, label=r'$\Theta_0^P(k=0.001/\text{Mpc})$')
        plt.plot(x, ThetaP_0_0001, label=r'$\Theta_0^P(k=0.0001/\text{Mpc})$')
        plt.plot(x, ThetaP_1_1, label=r'$\Theta_1^P(k=0.1/\text{Mpc})$')
        plt.plot(x, ThetaP_1_01, label=r'$\Theta_1^P(k=0.01/\text{Mpc})$')
        plt.plot(x, ThetaP_1_001, label=r'$\Theta_1^P(k=0.001/\text{Mpc})$')
        plt.plot(x, ThetaP_1_0001, label=r'$\Theta_1^P(k=0.0001/\text{Mpc})$')
        plt.plot(x, ThetaP_2_1, label=r'$\Theta_2^P(k=0.1/\text{Mpc})$')
        plt.plot(x, ThetaP_2_01, label=r'$\Theta_2^P(k=0.01/\text{Mpc})$')
        plt.plot(x, ThetaP_2_001, label=r'$\Theta_2^P(k=0.001/\text{Mpc})$')
        plt.plot(x, ThetaP_2_0001, label=r'$\Theta_2^P(k=0.0001/\text{Mpc})$')
        plt.xlabel(r'$x$')
        plt.ylabel(r'Polarization multipoles')
        plt.title(r'Evolution of polarization multipoles over time ($x$)')
        plt.legend()
        plt.grid(True)
        plt.savefig('../Plots/Milestone III/polarization.pdf')

def milestone3():
    
    
    
    #Run the functions.
    
    plot()