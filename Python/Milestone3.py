# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 12:06:38 2024

@author: Iván Villegas Pérez
"""

"""
This script contains all the relevant functions for analysing and plotting the
results of running the C++ codes related to solving the structure of the
Universe.

Here you can find the functions 'plot()' and
'milestone1()'. All of them are explained below.
"""

#Import all relevant packages

import numpy as np

import matplotlib.pyplot as plt

def plot():
    
    """
    This function generates and saves different plots based on the data
    obtained from 'perturbations_kX.txt' files, where X denotes the 'k=0.X'.
    The user is prompted to specify whether photon polarization and neutrinos 
    are included in the data for appropriate plotting. 
    The generated plots include density perturbations, velocity perturbations.
    and quadrupoles over time and, if polarization is included, polarizatoin
    multipoles.

    The function performs the following specific actions:
    1. Prompts the user to specify if photon polarization and neutrinos are
       included.
    2. Reads data from 'perturbations_kX.txt' files.
    3. Extracts relevant data such as time variable 'x', photon multipole
       values 'Theta_0_X', 'Theta_1_X', 'Theta_2_X', gravitational potential
       'Phi_X', and more.
    4. Generates plots for density perturbations, velocity perturbations, and 
       quadrupoles and, if polarization is included, polarizatoin
       multipoles using the extracted data.
    5. Saves the generated plots as PDF files in the '../Plots/Milestone III/'
       directory.

    Note: This function can be used to visualize the evolution of various
          cosmological perturbations over time and compare different components
          of the universe.

    Parameters:
        None.
    Returns:
        None.
    """
    
    #Check if photon polarization has been included.
    
    polarization_txt: str
    
    polarization_txt = input('Did you include photon polarization? [Yes, no]: ')
    
    polarization: bool
    
    if polarization_txt.lower()=='no':
        
        polarization = False
        
    else:
        
        polarization = True
        
    #Check if neutrinos have been included.
        
    neutrinos_txt: str
    
    neutrinos_txt = input('\nDid you include neutrinos? [Yes, no]: ')
    
    neutrinos: bool
    
    if neutrinos_txt.lower()=='no':
        
        neutrinos = False
        
    else:
        
        neutrinos = True
    
    #Read the data from 'perturbations_kX.txt' in folder 'Results',
    #where X denotes k=0.X.
    
    data_1 = np.loadtxt('../Results/perturbations_k1.txt')
    
    #x=ln(a): main time variable.
    
    x: np.array(float) = data_1[:, 0]
    
    #Theta_0_X: oreder 0 photon multipole.
    
    Theta_0_1: np.array(float) = data_1[:, 1]
    
    #Theta_1_X: oreder 1 photon multipole.
    
    Theta_1_1: np.array(float) = data_1[:, 2]
    
    #Theta_2_X: oreder 2 photon multipole.
    
    Theta_2_1: np.array(float) = data_1[:, 3]
    
    #Phi_X: gravitational potential.
    
    Phi_1: np.array(float) = data_1[:, 4]
    
    #Psi_X: curvature potential.
    
    Psi_1: np.array(float) = data_1[:, 5]
    
    #ThetaP_l_X: photon polarizaion multipoles.
    
    if polarization:
        ThetaP_0_1: np.array(float) = data_1[:, 6]
        ThetaP_1_1: np.array(float) = data_1[:, 7]
        ThetaP_2_1: np.array(float) = data_1[:, 8]
        
    #Nu_l_X: neutrino multipoles.
    
    if neutrinos:
        Nu_0_1: np.array(float) = data_1[:, 9]
        Nu_1_1: np.array(float) = data_1[:, 10]
        Nu_2_1: np.array(float) = data_1[:, 11]
        
    #delta_B_X: baryonic matter density perturbation
        
    delta_B_1: np.array(float) = data_1[:, 12]
    
    #delta_CDM_X: CDM density perturbation
        
    delta_CDM_1: np.array(float) = data_1[:, 13]
    
    #v_B_X: baryonic matter perturbation velocity
        
    v_B_1: np.array(float) = data_1[:, 14]
    
    #v_CDM_X: CDM perturbation velocity
        
    v_CDM_1: np.array(float) = data_1[:, 15]
    
    #New k
    
    data_01 = np.loadtxt('../Results/perturbations_k01.txt')
    Theta_0_01: np.array(float) = data_01[:, 1]
    Theta_1_01: np.array(float) = data_01[:, 2]
    Theta_2_01: np.array(float) = data_01[:, 3]
    Phi_01: np.array(float) = data_01[:, 4]
    Psi_01: np.array(float) = data_01[:, 5]
    if polarization:
        ThetaP_0_01: np.array(float) = data_01[:, 6]
        ThetaP_1_01: np.array(float) = data_01[:, 7]
        ThetaP_2_01: np.array(float) = data_01[:, 8]
    if neutrinos:
        Nu_0_01: np.array(float) = data_01[:, 9]
        Nu_1_01: np.array(float) = data_01[:, 10]
        Nu_2_01: np.array(float) = data_01[:, 11]
    delta_B_01: np.array(float) = data_01[:, 12]
    delta_CDM_01: np.array(float) = data_01[:, 13]
    v_B_01: np.array(float) = data_01[:, 14]
    v_CDM_01: np.array(float) = data_01[:, 15]
    
    #New k    
    
    data_001 = np.loadtxt('../Results/perturbations_k001.txt')    
    Theta_0_001: np.array(float) = data_001[:, 1]
    Theta_1_001: np.array(float) = data_001[:, 2]
    Theta_2_001: np.array(float) = data_001[:, 3]
    Phi_001: np.array(float) = data_001[:, 4]
    Psi_001: np.array(float) = data_001[:, 5]
    if polarization:
        ThetaP_0_001: np.array(float) = data_001[:, 6]
        ThetaP_1_001: np.array(float) = data_001[:, 7]
        ThetaP_2_001: np.array(float) = data_001[:, 8]
    if neutrinos:
        Nu_0_001: np.array(float) = data_001[:, 9]
        Nu_1_001: np.array(float) = data_001[:, 10]
        Nu_2_001: np.array(float) = data_001[:, 11]        
    delta_B_001: np.array(float) = data_001[:, 12]
    delta_CDM_001: np.array(float) = data_001[:, 13]
    v_B_001: np.array(float) = data_001[:, 14]
    v_CDM_001: np.array(float) = data_001[:, 15]
    
    #New k    
        
    data_0001 = np.loadtxt('../Results/perturbations_k0001.txt')    
    Theta_0_0001: np.array(float) = data_0001[:, 1]
    Theta_1_0001: np.array(float) = data_0001[:, 2]
    Theta_2_0001: np.array(float) = data_0001[:, 3]
    Phi_0001: np.array(float) = data_0001[:, 4]
    Psi_0001: np.array(float) = data_0001[:, 5]
    if polarization:
        ThetaP_0_0001: np.array(float) = data_0001[:, 6]
        ThetaP_1_0001: np.array(float) = data_0001[:, 7]
        ThetaP_2_0001: np.array(float) = data_0001[:, 8]
    if neutrinos:
        Nu_0_0001: np.array(float) = data_0001[:, 9]
        Nu_1_0001: np.array(float) = data_0001[:, 10]
        Nu_2_0001: np.array(float) = data_0001[:, 11]
    delta_B_0001: np.array(float) = data_0001[:, 12]
    delta_CDM_0001: np.array(float) = data_0001[:, 13]
    v_B_0001: np.array(float) = data_0001[:, 14]
    v_CDM_0001: np.array(float) = data_0001[:, 15]
    
    
    #Make the different plots
    
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
    plt.yscale('log')
    plt.xlabel(r'$x$')
    plt.ylabel(r'Density perturbation ($\delta_i$)')
    plt.title(r'Evolution of density perturbations over time ($x$)')
    plt.legend()
    plt.grid(True)
    plt.xlim(x[0], x[-1])
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
    plt.yscale('log')
    plt.xlabel(r'$x$')
    plt.ylabel(r'Velocity perturbation ($v_i$)')
    plt.title(r'Evolution of velocity perturbations over time ($x$)')
    plt.legend()
    plt.grid(True)
    plt.xlim(x[0], x[-1])
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
    plt.xlim(x[0], x[-1])
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
    plt.xlim(x[0], x[-1])
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
        plt.xlim(x[0], x[-1])
        plt.savefig('../Plots/Milestone III/polarization.pdf')

def milestone3():
    
    """
    This function serves as the entry point for running the tasks and functions
    related to Milestone III of the cosmology research project.
    In this milestone, the 'plot' function is called to generate and save plots
    based on the data obtained from 'perturbations_kX.txt' file.

    The function performs the following specific actions:
    1. Calls the 'plot' function to generate and save plots for density
    perturbations, velocity perturbations quadrupoles and, if polarization is 
    included, polarization multipoles over time.
    
    Note: The purpose of this function is to execute the necessary tasks and
    functions for Milestone III of the cosmology research project. It can be
    customized to include additional steps or functions relevant to this 
    milestone.

    Parameters:
        None.
    Returns:
        None.
    """
    
    #Run the functions.
    
    plot()
