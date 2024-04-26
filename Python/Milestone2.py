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

def plots(neutrinos: bool):
    
    """
    Generates plots from the data in 'recombination.txt' and saves them as PDF
    files.

    Parameters:
        None.

    Returns:
        None.

    Reads the data from 'recombination.txt' or 'recombination.txt' (if
    neutrinos are included) in the 'Results' folder, assigns the data to
    corresponding variables, and creates various plots from the data, saving
    them as PDF files.  It calls the 'table()' function to print a table of
    relevant data.
    """
    
    #Read the data from 'recombination.txt' in folder 'Results'.
    
    if neutrinos:
        data = np.loadtxt('../Results/recombination_neutrinos.txt')
    else:
        data = np.loadtxt('../Results/recombination.txt')
    
    #Assign the different variables to its corresponding data.
    
    #x=ln(a): main time variable.
    
    x: np.array(float) = data[:, 0]
    
    #X_e: fractional electron density.
    
    X_e: np.array(float) = data[:, 1]
    
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
    
    #Indeces for decoupling and last scattering
    
    last_scattering: int
    
    recombination: int
    
    last_scattering, recombination = tables(neutrinos)
    
    #Make the different plots and save them.
    
    plt.figure()
    plt.plot(x, X_e, label=r'$X_e$')
    plt.plot(x, X_e_Saha, label=r'$X_e$ from Saha equation', ls='dashed')
    plt.hlines(X_e[-1]-1e-5, x[0], x[-1], ls='dotted',\
               label=r'$X_e$(Freeze-out)$\approx$'+f'{X_e[-1]*1e4:.3f}'+\
                   r'$\cdot10^{-4}$', color='red')
    plt.vlines(x[recombination], 1e-4, 2, ls='dotted', label='Recombination',\
               color='purple')
    plt.grid(True)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$X_e$')
    plt.title(r'$X_e$ against $x$')
    plt.yscale("log")
    plt.xlim(x[0], x[-1])
    plt.ylim(1e-4, 2)
    plt.legend()
    if neutrinos:
        plt.savefig('../Plots/Milestone II/X_e_neutrinos.pdf')
    else:
        plt.savefig('../Plots/Milestone II/X_e.pdf')
    
    plt.figure()
    plt.plot(x, tau, label=r'$\tau$')
    plt.plot(x, -dtau, label=r'$-\frac{\text{d}\tau}{\text{d}x}$')
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
    plt.plot(x, dg_tilde/max(abs(dg_tilde)),\
             label=r'$\frac{\text{d}\tilde{g}}{\text{d}x}$'+\
                   r'/|$\frac{\text{d}\tilde{g}}{\text{d}x}$|$_\text{max}$')
    plt.plot(x, ddg_tilde/max(abs(ddg_tilde)),\
             label=r'$\frac{\text{d}^2\tilde{g}}{\text{d}x^2}$'+\
                   r'/|$\frac{\text{d}^2\tilde{g}}{\text{d}x^2}$|$_\text{max}$')
    plt.vlines(x[last_scattering], -1.1, 1.1, ls='dashed',\
               label='Last scattering surface', color='black')
    plt.xlabel(r'x')
    plt.title('Visibility function and its evolution')
    plt.grid(True)
    plt.legend()
    plt.ylim(-1.1, 1.1)
    plt.xlim(-7.4, -6.4)
    plt.savefig('../Plots/Milestone II/g_tilde and derivatives.pdf')
    
def tables(neutrinos: bool):
    
    """
    Generates a table from the data in 'recombination.txt' or
    'recombination.txt' (if neutrinos are included) and prints it.

    Parameters:
        None.

    Returns:
        None.

    Reads the data from 'recombination.txt' or 'recombination.txt' (if
    neutrinos are included) in the 'Results' folder, assigns the data to
    corresponding variables, looks for specific value in the data and
    creates/prints a table summarizing these data.
    """
    
    #Read the data from 'cosmology.txt' in folder 'Results'.
    
    if neutrinos:
        data = np.loadtxt('../Results/recombination_neutrinos.txt')
    else:
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
    
    #t: cosmological time.
    
    t: np.array(float) = data[:, 10]/(1e9*365*24*60*60)
    
    #z: redshift.
    
    z: np.array(float) = np.exp(-x)-1
    
    #Find indices for decoupling, last scattering, half-way recombination,
    #half-way recombination using the Saha approximation and and end of the
    #recombination time.
    
    decoupling: int = int(input('\nIndex at decoupling: '))
    
    last_scattering: int = np.argmin(np.abs(tau-1))
    
    half_way: int = np.argmin(np.abs(X_e-0.5))
    
    half_way_Saha: int = np.argmin(np.abs(X_e_Saha-0.5))
    
    recombination: int = np.argmin(np.abs(X_e-0.1))
    
    print(f'\nTIMES\n\
                            x       z        t [Myr]\n\
    Decoupling:           {x[decoupling]:.3f}  {z[decoupling]:.2f}     {t[decoupling]*1e3:.3f}\n\
    Last scattering:      {x[last_scattering]:.3f}  {z[last_scattering]:.2f}     {t[last_scattering]*1e3:.3f}\n\
    Half-way rec:         {x[half_way]:.3f}  {z[half_way]:.2f}     {t[half_way]*1e3:.3f}\n\
    Half-way rec (Saha):  {x[half_way_Saha]:.3f}  {z[half_way_Saha]:.2f}     {t[half_way_Saha]*1e3:.3f}\n\
    Recombination:        {x[recombination]:.3f}  {z[recombination]:.2f}     {t[recombination]*1e3:.3f}')
          
    return(last_scattering, recombination)
    
def milestone2(neutrinos: bool):
    
    """
    Execute main functionality for recombination history solving.

    This function serves as the entry point for the recombination history
    program. It calls the 'plots()' function to generate plots.
    
    Parameters:
        None.

    Returns:
        None.
    """
    
    #Run the functions.
    
    plots(neutrinos)
    