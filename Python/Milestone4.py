# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 12:06:38 2024

@author: Iván Villegas Pérez
"""

#Import all relevant packages

import scipy as sc

import numpy as np

import matplotlib.pyplot as plt

def plot(polarization: bool):
    
    data = np.loadtxt('../Results/cells.txt')
    
    ell: np.array(float) = data[:, 0]
    
    C_ell: np.array(float) = data[:, 1]
    
    data_ACT = np.loadtxt('../Data/ACT_data.txt', skiprows=1)
    
    k_ACT: np.array(float) = data_ACT[:, 0]
    
    P_ACT: np.array(float) = data_ACT[:, 1]
    
    P_upper: np.array(float) = data_ACT[:, 2]
    
    data_CMB_PS = np.loadtxt('../Data/CMB_PS_data.txt', skiprows=1)
    
    k_CMB_PS: np.array(float) = data_CMB_PS[:, 0]
    
    P_CMB_PS: np.array(float) = data_CMB_PS[:, 1]
    
    ErrorP_CMB_PS: np.array(float) = data_CMB_PS[:, 2]
    
    if polarization:
    
        data_highEE = np.loadtxt('../Data/high_EE_data.txt', skiprows=1)
        
        l_highEE: np.array(float) = data_highEE[:, 0]
        
        D_l_highEE: np.array(float) = data_highEE[:, 1]
        
        DeltaD_down_highEE: np.array(float) = data_highEE[:, 2]
        
        DeltaD_up_highEE: np.array(float) = data_highEE[:, 3]
        
        data_highTE = np.loadtxt('../Data/high_TE_data.txt', skiprows=1)
        
        l_highTE: np.array(float) = data_highTE[:, 0]
        
        D_l_highTE: np.array(float) = data_highTE[:, 1]
        
        DeltaD_down_highTE: np.array(float) = data_highTE[:, 2]
        
        DeltaD_up_highTE: np.array(float) = data_highTE[:, 3]
        
        data_highTT = np.loadtxt('../Data/high_TT_data.txt', skiprows=1)
        
        l_highTT: np.array(float) = data_highTT[:, 0]
        
        D_l_highTT: np.array(float) = data_highTT[:, 1]
        
        DeltaD_down_highTT: np.array(float) = data_highTT[:, 2]
        
        DeltaD_up_highTT: np.array(float) = data_highTT[:, 3]
    
    data_lowTT = np.loadtxt('../Data/low_TT_data.txt', skiprows=1)
    
    l_lowTT: np.array(float) = data_lowTT[:, 0]
    
    D_l_lowTT: np.array(float) = data_lowTT[:, 1]
    
    DeltaD_down_lowTT: np.array(float) = data_lowTT[:, 2]
    
    DeltaD_up_lowTT: np.array(float) = data_lowTT[:, 3]
    
    
def milestone4(polarization: bool):
    
    plot(polarization)
    