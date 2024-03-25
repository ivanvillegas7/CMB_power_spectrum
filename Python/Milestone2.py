# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 12:05:37 2024

@author: Iván Villegas Pérez
"""

#Import all relevant packages

import scipy as sc

import numpy as np

import matplotlib.pyplot as plt

def milestone2():
    
    #Read the data from 'cosmology.txt' in folder 'Results'.
    
    data = np.loadtxt('../Results/recombination.txt')
    
    #Assign the different variables to its corresponding data.
    
    #x=ln(a): main time variable.
    
    x: np.array(float) = data[:, 0]
    
    #X_e: fractional electron density.
    
    X_e: np.array(float) = data[:, 1]
    
    #n_e: electron density number
    
    n_e: np.array(float) = data[:, 2]
    
    tau: np.array(float) = data[:, 3]
    
    #dtau: first derivative of tau with respect to x.
    
    dtau: np.array(float) = data[:, 4]
    
    #ddtau: second derivative of tau with respect to x.
    
    ddtau: np.array(float) = data[:, 5]
    
    g_tilde: np.array(float) = data[:, 6]
    
    #dg_tilde: first derivative of g_tilde with respect to x.
    
    dg_tilde: np.array(float) = data[:, 7]
    
    #ddg_tilde: second derivative of g_tilde with respect to x.
    
    ddg_tilde: np.array(float) = data[:, 8]
    
    #Make the different plots and save them.
    
    plt.figure()
    plt.plot(x, X_e)
    plt.grid(True)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$X_e$')
    plt.title(r'$X_e$ against $x$')
    plt.savefig('../Plots/Milestone II/X_e.pdf')
    
    plt.figure()
    plt.plot(x, tau, label=r'$\tau$')
    plt.plot(x, dtau, label=r'$\frac{\text{d}\tau}{\text{d}x}$')
    plt.plot(x, ddtau, label=r'$\frac{\text{d}^2\tau}{\text{d}x^2}$')
    plt.xlabel(r'x')
    plt.ylabel(r'$\tau$, $\frac{\text{d}\tau}{\text{d}x}$ and $\frac{\text{d}^2\tau}{\text{d}x^2}$')
    plt.title(r'$\tau$, $\frac{\text{d}\tau}{\text{d}x}$ and $\frac{\text{d}^2\tau}{\text{d}x^2}$ against $x$')
    plt.yscale("log")
    plt.grid(True)
    plt.legend()
    plt.savefig('../Plots/Milestone II/tau and derivatives.pdf')
    
    plt.figure()
    plt.plot(x, g_tilde, label=r'$\tilde{g}$')
    plt.plot(x, dg_tilde, label=r'$\frac{\text{d}\tilde{g}}{\text{d}x}$')
    plt.plot(x, ddg_tilde, label=r'$\frac{\text{d}^2\tilde{g}}{\text{d}x^2}$')
    plt.xlabel(r'x')
    plt.ylabel(r'$\tilde{g}$, $\frac{\text{d}\tilde{g}}{\text{d}x}$ and $\frac{\text{d}^2\tilde{g}}{\text{d}x^2}$')
    plt.title(r'$\tilde{g}$, $\frac{\text{d}\tilde{g}}{\text{d}x}$ and $\frac{\text{d}^2\tilde{g}}{\text{d}x^2}$ against $x$')
    plt.grid(True)
    plt.legend()
    plt.savefig('../Plots/Milestone II/g_tilde and derivatives.pdf')
    