#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 11:20:46 2024

@author: Iván Villegas Pérez
"""

"""
This script contains all the auxiliar functions for analysing and plotting the
results of running the C++ codes related to computing the CMB power spectrum.

Here you can find the functions 'plot()' and 'E()', 'D_L()' and 'find_value()'.
All of them are explained below.
"""

import scipy as sc

import numpy as np

import os

import matplotlib

#Whether to also display figures in a floating window, in addition to saving
#them to disk as always. OFF by default so that running './cmb' stays fully
#automatic (no windows popping up, nothing to close by hand). Turn it on
#with:
#
#   CMB_SHOW_PLOTS=1 ./cmb
#
#or, to only view the plots of a single milestone without rerunning the C++
#part:
#
#   CMB_SHOW_PLOTS=1 python3 Milestone4.py

SHOW_PLOTS: bool = os.environ.get('CMB_SHOW_PLOTS', '0') == '1'

#Pick the backend BEFORE 'matplotlib.pyplot' gets imported anywhere (the
#backend can only be chosen once per process). When no window is going to be
#shown, 'Agg' is used: a purely headless, file-only backend that never talks
#to Qt/GTK/Wayland, so no GUI warning of any kind can be printed. When a
#window IS wanted, 'TkAgg' (Tkinter, bundled with Python, no Qt involved) is
#used instead of the default Qt backend, which is what was printing the
#"Ignoring XDG_SESSION_TYPE=wayland" warning on Gnome/Wayland.

matplotlib.use('Agg' if not SHOW_PLOTS else 'TkAgg')

import matplotlib.pyplot as plt

from typing import List

def show_or_save(filename: str):

    """
    Save the current matplotlib figure to 'filename', exactly like
    plt.savefig() always did. If the environment variable CMB_SHOW_PLOTS=1
    is set, the figure is also kept open to be displayed in a floating
    window later (see finalize_plots()), instead of being closed straight
    away.

    Parameters:
        filename : str
            Path (including extension) to save the figure to.

    Returns:
        None.
    """

    plt.savefig(filename, bbox_inches='tight')

    if not SHOW_PLOTS:

        #Free the figure's memory immediately, same behaviour as before.

        plt.close()

def finalize_plots():

    """
    Show every figure created since the last call to finalize_plots(), each
    in its own floating window, and block until all of them have been
    closed by hand. Does nothing (and closes nothing) if CMB_SHOW_PLOTS is
    not set to '1', so the default automated pipeline is unaffected.

    Call this once, at the very end of each milestoneX() function, after
    all of that milestone's plots have already been created/saved.

    Parameters:
        None.

    Returns:
        None.
    """

    if SHOW_PLOTS and plt.get_fignums():

        plt.show()

#Define function 'plot()'.

def plot(x: np.array(float), y: np.array(float), i1: int, i2:int):
    
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
    
def E(z: float, H0: float, Omega_M: float, Omega_K: float, Omega_R: float,\
      Omega_Lambda: float):

    return np.sqrt(Omega_M*(1+z)**3+Omega_K*(1+z)**2+Omega_R*(1+z)**4+\
                   Omega_Lambda)
    
def D_L(z: np.array(float), H0: float, Omega_M: float, Omega_K: float,\
        Omega_R: float, Omega_Lambda: float):

    """
    Calculates the luminosity distance for a given redshift using the
    cosmological parameters.

    Parameters:
        z : np.array(float)
            Array of redshift values.
        H0 : float
            Hubble constant at present time.
        Omega_M : float
            Matter density parameter.
        Omega_K : float
            Curvature density parameter.
        Omega_R : float
            Radiation density parameter.
        Omega_Lambda : float
            Dark energy density parameter.

    Returns:
        d_L : List[float]
            Luminosity distance corresponding to the given redshifts.

    The function calculates the luminosity distance using the given redshift
    values and the cosmological parameters. It computes the E(z) term and then
    integrates 1/E(z) to find the corresponding luminosity distance for each
    redshift value. The resulting luminosity distances are returned as a list.
    """
    
    d_L: List[float] = []
    
    for i in range(len(z)):

        x: np.array = np.array(np.linspace(0, z[i]))
    
        d_L.append(((sc.constants.c/H0)*(1+z[i])\
                    *sc.integrate.simpson(1/E(x, H0, Omega_M, Omega_K, Omega_R,\
                                            Omega_Lambda), x)))
        
    return np.array(d_L)

def read_flags(path: str = '../Results/run_info.txt'):

    """
    Read the flags (neutrinos, Helium, reionization, polarization) and
    derived scalars (decoupling index, sound horizon, ...) written by the
    C++ code to 'Results/run_info.txt'.

    This avoids having to copy values printed to the terminal by the C++
    code (e.g. 'Index at decoupling') by hand into the Python scripts.

    Parameters:
        path : str
            Path to the file written by RecombinationHistory::output_info().

    Returns:
        dict
            Dictionary mapping each flag/value name to an int or float.
    """

    flags = {}

    with open(path) as f:

        for line in f:

            line = line.strip()

            if not line:

                continue

            key, value = line.split(maxsplit=1)

            try:

                flags[key] = int(value)

            except ValueError:

                flags[key] = float(value)

    return flags

def index_equality():
    
    """
    Find indices where the densities of cold matter and relativistic particles
    are approximately equal.

    This function reads data from a file containing cosmological information
    and calculates the indices where the density of cold matter becomes similar
    to the density of relativistic particles.
    
    Parameters:
        None.

    Returns:
        Tuple of integers (index_M_R, index_M_Lambda).
    """
    
    #Read the data from 'cosmology.txt' in folder 'Results'.
    
    data = np.loadtxt('../Results/cosmology.txt')
    
    #OmegaLambda: relative density of dark energy.
    
    OmegaLambda: np.array(float) = data[:, 6]
    
    #Define the omegas corresponding to cold and relativistic matter/particles.
        
    OmegaM:np.array(float) = data[:, 4]+data[:, 5]

    OmegaRel:np.array(float) = data[:, 7]+data[:, 8]
    
    #Get the indeces when the cuantity of cold matter is the same (or most 
    #similar) to the cuantity of relativistic particles.
    
    index_M_R: int = np.argmin(np.abs(OmegaRel-OmegaM))
    
    index_M_Lambda: int = index_M_R+\
                          np.argmin(np.abs(OmegaLambda[index_M_R:]\
                                           -OmegaM[index_M_R:]))
                              
    return(index_M_R, index_M_Lambda)