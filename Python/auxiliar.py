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

import matplotlib.pyplot as plt

from typing import List

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
    integrates to find the corresponding luminosity distance for each redshift
    value. The resulting luminosity distances are returned as a list.
    """
    
    d_L: List[float] = []
    
    for i in range(len(z)):

        x: np.array = np.array(np.linspace(0, z[i]))
    
        d_L.append(((sc.constants.c/H0)*(1+z[i])\
                    *sc.integrate.simpson(E(x, H0, Omega_M, Omega_K, Omega_R,\
                                            Omega_Lambda), x)))
        
    return np.array(d_L)
