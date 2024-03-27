# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 12:05:37 2024

@author: Iván Villegas Pérez
"""

"""
This script contains all the relevant functions for analysing and plotting the
results of running the C++ codes related to solving the recombination history
of the Universe.

Here you can find the functions 'plots()' and 'milestone2()'.
All of them are explained below.
"""

#Import all relevant packages

import numpy as np

import matplotlib.pyplot as plt

def plots():
    
    """
    Generates plots from the data in 'recombination.txt' and saves them as PDF
    files.

    Parameters:
        None.

    Returns:
        None.

    Reads the data from 'recombination.txt' in the 'Results' folder, assigns
    the data to corresponding variables, and creates various plots from the
    data, saving them as PDF files.
    """
    
    #Read the data from 'cosmology.txt' in folder 'Results'.
    
    data = np.loadtxt('../Results/recombination.txt')
    
    #Assign the different variables to its corresponding data.
    
    #x=ln(a): main time variable.
    
    x: np.array(float) = data[:, 0]
    
    #X_e: fractional electron density.
    
    X_e: np.array(float) = data[:, 1]
    
    #n_e: electron density number
    
    n_e: np.array(float) = data[:, 2]
    
    #tau: optical depth
    
    tau: np.array(float) = data[:, 3]
    
    #dtau: first derivative of tau with respect to x.
    
    dtau: np.array(float) = data[:, 4]
    
    #ddtau: second derivative of tau with respect to x.
    
    ddtau: np.array(float) = data[:, 5]
    
    #g_tilde: visibility function
    
    g_tilde: np.array(float) = data[:, 6]
    
    #dg_tilde: first derivative of g_tilde with respect to x.
    
    dg_tilde: np.array(float) = data[:, 7]
    
    #ddg_tilde: second derivative of g_tilde with respect to x.
    
    ddg_tilde: np.array(float) = data[:, 8]
    
    #X_e_Saha: fractional electron density only from Saha equation.
    
    X_e_Saha: np.array(float) = data[:, 9]
    
    #Make the different plots and save them.
    
    plt.figure()
    plt.plot(x, X_e, label=r'$X_e$')
    plt.plot(x, X_e_Saha, label=r'$X_e$ from Saha equation', ls='dashed')
    plt.grid(True)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$X_e$')
    plt.title(r'$X_e$ against $x$')
    plt.yscale("log")
    plt.ylim(1e-4, 2)
    plt.savefig('../Plots/Milestone II/X_e.pdf')
    
    plt.figure()
    plt.plot(x, tau, label=r'$\tau$')
    plt.plot(x, dtau, label=r'$\frac{\text{d}\tau}{\text{d}x}$')
    plt.plot(x, ddtau, label=r'$\frac{\text{d}^2\tau}{\text{d}x^2}$')
    plt.xlabel(r'x')
    plt.title('Optical depth and its evolution')
    plt.yscale("log")
    plt.grid(True)
    plt.legend()
    plt.savefig('../Plots/Milestone II/tau and derivatives.pdf')
    
    plt.figure()
    plt.plot(x, g_tilde/max(abs(g_tilde)), label=r'$\tilde{g}$'+\
                                                 r'/|$\tilde{g}$|$_\text{max}$')
    plt.plot(x, dg_tilde/max(abs(dg_tilde)), ls='dashed',\
             label=r'$\frac{\text{d}\tilde{g}}{\text{d}x}$'+\
                   r'/|$\frac{\text{d}\tilde{g}}{\text{d}x}$|$_\text{max}$')
    plt.plot(x, ddg_tilde/max(abs(ddg_tilde)), ls='dotted',\
             label=r'$\frac{\text{d}^2\tilde{g}}{\text{d}x^2}$'+\
                   r'/|$\frac{\text{d}^2\tilde{g}}{\text{d}x^2}$|$_\text{max}$')
    plt.xlabel(r'x')
    plt.title('Visibility function and its evolution')
    plt.grid(True)
    plt.legend()
    plt.xlim(-6.94, -6.88)
    plt.savefig('../Plots/Milestone II/g_tilde and derivatives.pdf')
    
def tables():
    
    """
    Generates a table from the data in 'recombination.txt' and prints it.

    Parameters:
        None.

    Returns:
        None.

    Reads the data from 'recombination.txt' in the 'Results' folder, assigns
    the data to corresponding variables, looks for specific value in the data
    and creates/prints a table summarizing these data.
    """
    
    #Read the data from 'cosmology.txt' in folder 'Results'.
    
    data = np.loadtxt('../Results/recombination.txt')
    
    #Assign the different variables to its corresponding data.
    
    #x=ln(a): main time variable.
    
    x: np.array(float) = data[:, 0]
    
    #X_e: fractional electron density.
    
    X_e: np.array(float) = data[:, 1]
    
    #tau: optical depth
    
    tau: np.array(float) = data[:, 3]
    
    #X_e_Saha: fractional electron density only from Saha equation.
    
    X_e_Saha: np.array(float) = data[:, 9]
    
    
    
def milestone2():
    
    """
    Execute main functionality for recombination history solving.

    This function serves as the entry point for the recombination history
    program. It calls the 'plots()' function to generate plots and the
    'table()' function to print a table of relevant data.
    
    Parameters:
        None.

    Returns:
        None.
    """
    
    #Run the functions.
    
    plots()
    
    tables()
    