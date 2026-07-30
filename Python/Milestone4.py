# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 12:06:38 2024

@author: Iván Villegas Pérez
"""

"""
This script contains all the relevant functions for plotting the results of
running the C++ codes related to computing the CMB power-spectrum.

Here you can find the functions 'CMB_PS_check()', 'CMB_PowerSpectrum()',
'Matter_PowerSpectrum()' 'other_plots()', 'CMB_map()' and 'milestone4()'.
All of them are explained below.
"""

#Import all relevant packages

import scipy as sc

import numpy as np

import io

import auxiliar as aux

# auxiliar picks the matplotlib backend (Agg / TkAgg) before pyplot is
# imported anywhere, so import it first to avoid GUI backend warnings.
import matplotlib.pyplot as plt

try:
    import notification as notif
except ImportError:
    notif = None

#JPG snapshot of the CMB TT power-spectrum figure, captured in
#CMB_PowerSpectrum() and sent as part of the "simulation finished"
#Telegram notification at the end of milestone4().

_last_cmb_ps_jpg = None
_last_cmb_ps_filename = None

def CMB_PS_components():
    
    """
    Check and plot various components of the CMB power-spectrum.

    This function reads data from different files containing the CMB
    power-spectrum information for different components such as SW, ISW,
    Doppler, and Polarization. It then generates a plot to compare these
    components with the overall theoretical prediction of the CMB
    power-spectrum.

    Parameters:
        None.

    Returns:
        None.
    """
    
    #Read the data from 'cells.txt' in folder 'Results'
    
    data_Cells = np.loadtxt('../Results/cells.txt')
    
    #ell: multipole moment.
    
    ell: np.array(float) = data_Cells[:, 0]
    
    #C_ell: power-spectrum.
    
    C_ell_TT: np.array(float) = data_Cells[:, 1]
    
    #Read the data from 'cells_SW.txt' in folder 'Results'
    
    data_Cells_SW = np.loadtxt('../Results/cells_SW.txt')
    
    #C_ell: power-spectrum.
    
    C_ell_TT_SW: np.array(float) = data_Cells_SW[:, 1]
    
    #Read the data from 'cells_ISW.txt' in folder 'Results'
    
    data_Cells_ISW = np.loadtxt('../Results/cells_ISW.txt')
    
    #C_ell: power-spectrum.
    
    C_ell_TT_ISW: np.array(float) = data_Cells_ISW[:, 1]
    
    #Read the data from 'cells_DOPPLER.txt' in folder 'Results'
    
    data_Cells_DOPPLER = np.loadtxt('../Results/cells_DOPPLER.txt')
    
    #C_ell: power-spectrum.
    
    C_ell_TT_DOPPLER: np.array(float) = data_Cells_DOPPLER[:, 1]
    
    #Read the data from 'cells_POL.txt' in folder 'Results'
    
    data_Cells_POL = np.loadtxt('../Results/cells_POL.txt')
    
    #C_ell: power-spectrum.
    
    C_ell_TT_POL: np.array(float) = data_Cells_POL[:, 1]
    
    #Make the plot.
        
    plt.figure()
    plt.plot(ell, C_ell_TT, label='Theory prediction', ls='solid')
    plt.plot(ell, C_ell_TT_SW, label='SW theory prediction', ls='dashed')
    plt.plot(ell, C_ell_TT_ISW, label='ISW theory prediction', ls='dashed')
    plt.plot(ell, C_ell_TT_DOPPLER, label='Doppler theory prediction',\
             ls='dashed')
    plt.plot(ell, C_ell_TT_POL, label='Polarization theory prediction',\
             ls='dashed')
    plt.xlabel(r'Multipole $\ell$')
    plt.ylabel(r'$\ell(\ell+1)C_\ell/2\pi$ [$\mu$K$^2$]')
    plt.title('CMB power-spectrum')
    plt.xscale('log')
    plt.legend()
    plt.grid()
    aux.show_or_save('../Plots/Milestone IV/CMB PS components.pdf')

def CMB_PowerSpectrum():
    
    """
    Plot the Cosmic Microwave Background (CMB) power-spectrum.

    This function reads data from a file containing CMB power-spectrum
    information, including multipole moments (ell) and their corresponding
    power values (C_ell). It then generates a plot to visualize the theoretical
    prediction of the CMB power-spectrum, along with observational data.
    
    Parameters:
        None.

    Returns:
        None.
    """
    
    #Check if polarization has been included, reading the flag written by
    #the C++ code instead of asking the user.
    
    polarization: bool = bool(aux.read_flags()['polarization'])
    
    #Read the data from 'cells.txt' and 'Matter_PS.txt' in folder 'Results'
    
    data_Cells = np.loadtxt('../Results/cells.txt')
    
    #ell: multipole moment.
    
    ell: np.array(float) = data_Cells[:, 0]
    
    #C_ell: power-spectrum.
    
    C_ell_TT: np.array(float) = data_Cells[:, 1]
    
    #Read and assign the data from 'low_TT.txt' in folder 'Data'.
    
    data_lowTT = np.loadtxt('../Data/low_TT.txt', skiprows=1)
    
    l_lowTT: np.array(float) = data_lowTT[:, 0]
    
    D_l_lowTT: np.array(float) = data_lowTT[:, 1]
    
    DeltaD_down_lowTT: np.array(float) = data_lowTT[:, 2]
    
    DeltaD_up_lowTT: np.array(float) = data_lowTT[:, 3]
    
    #Read and assign the data from 'high_TT.txt' in folder 'Data'.
    
    data_highTT = np.loadtxt('../Data/high_TT.txt', skiprows=1)
    
    l_highTT: np.array(float) = data_highTT[:, 0]
    
    D_l_highTT: np.array(float) = data_highTT[:, 1]
    
    DeltaD_down_highTT: np.array(float) = data_highTT[:, 2]
    
    DeltaD_up_highTT: np.array(float) = data_highTT[:, 3]
    
    if polarization:

        C_ell_EE: np.array(float) = data_Cells[:, 2]

        C_ell_TE: np.array(float) = data_Cells[:, 3]
        
        #Read and assign the data from 'high_EE.txt' in folder 'Data'.
    
        data_highEE = np.loadtxt('../Data/high_EE.txt', skiprows=1)
        
        l_highEE: np.array(float) = data_highEE[:, 0]
        
        D_l_highEE: np.array(float) = data_highEE[:, 1]
        
        DeltaD_down_highEE: np.array(float) = data_highEE[:, 2]
        
        DeltaD_up_highEE: np.array(float) = data_highEE[:, 3]
        
        #Read and assign the data from 'high_TE.txt' in folder 'Data'.
        
        data_highTE = np.loadtxt('../Data/high_TE.txt', skiprows=1)
        
        l_highTE: np.array(float) = data_highTE[:, 0]
        
        D_l_highTE: np.array(float) = data_highTE[:, 1]
        
        DeltaD_down_highTE: np.array(float) = data_highTE[:, 2]
        
        DeltaD_up_highTE: np.array(float) = data_highTE[:, 3]
             
    #Make the plot.
        
    plt.figure()
    plt.plot(ell, C_ell_TT, label='Theory prediction', ls='solid')
    plt.errorbar(l_lowTT, D_l_lowTT, yerr=[DeltaD_down_lowTT, DeltaD_up_lowTT],\
                ls='none', label=r'Low $\ell$ TT data', marker='.', capsize=2)
    plt.errorbar(l_highTT, D_l_highTT,\
                yerr=[DeltaD_down_highTT, DeltaD_up_highTT], ls='none',\
                label=r'High $\ell$ TT data', marker='.', capsize=2)
    plt.xlabel(r'Multipole $\ell$')
    plt.ylabel(r'$\ell(\ell+1)C_\ell/2\pi$ [$\mu$K$^2$]')
    plt.title('CMB power-spectrum')
    plt.xscale('log')
    plt.legend()
    plt.grid()
    
    #Capture a JPG snapshot of this figure for the Telegram notification
    #sent at the end of the whole run (see milestone4()), before
    #aux.show_or_save() potentially closes it.
    
    global _last_cmb_ps_jpg, _last_cmb_ps_filename
    
    _cmb_ps_buffer = io.BytesIO()
    
    try:
        
        #Saving as JPG requires Pillow; fall back to PNG if it's missing.
        
        plt.savefig(_cmb_ps_buffer, format='jpg', bbox_inches='tight', facecolor='white')
        
        _last_cmb_ps_filename = 'CMB_power_spectrum.jpg'
        
    except Exception:
        
        _cmb_ps_buffer = io.BytesIO()
        
        plt.savefig(_cmb_ps_buffer, format='png', bbox_inches='tight', facecolor='white')
        
        _last_cmb_ps_filename = 'CMB_power_spectrum.png'
    
    _cmb_ps_buffer.seek(0)
    
    _last_cmb_ps_jpg = _cmb_ps_buffer.getvalue()
    
    aux.show_or_save('../Plots/Milestone IV/CMB PS.pdf')
    
    if polarization:
        
        plt.figure()
        plt.plot(ell, C_ell_TE, label='Theory prediction', ls='solid')
        plt.errorbar(l_highTE, D_l_highTE,\
                     yerr=[DeltaD_down_highTE, DeltaD_up_highTE], ls='none',\
                     label=r'High $\ell$ TE data', marker='.', capsize=2)
        plt.xlabel(r'Multipole $\ell$')
        plt.ylabel(r'$\ell(\ell+1)C_\ell^{TE}/2\pi$ [$\mu$K$^2$]')
        plt.title('The temperature-polarization cross power-spectrum')
        plt.legend()
        plt.grid()
        aux.show_or_save('../Plots/Milestone IV/CMB_TE PS.pdf')
        
        plt.figure()
        plt.plot(ell, C_ell_EE, label='Theory prediction')
        plt.errorbar(l_highEE, D_l_highEE,\
                     yerr=[DeltaD_down_highEE, DeltaD_up_highEE], ls='none',\
                     label=r'High $\ell$ EE data', marker='.', capsize=2)
        plt.xlabel(r'Multipole $\ell$')
        plt.ylabel(r'$\ell(\ell+1)C_\ell^{TE}/2\pi$ [$\mu$K$^2$]')
        plt.title('The (E mode) polarization power-spectrum')
        plt.legend()
        plt.grid()
        aux.show_or_save('../Plots/Milestone IV/CMB_EE PS.pdf')
    
def Matter_PowerSpectrum():
    
    """
    Plot the matter power spectrum.

    This function reads data from files containing matter power-spectrum
    information,as well as data from various surveys and experiments.
    It then generates a plot to visualize the theoretical prediction of the
    matter power-spectrum along with observational data.

    Parameters:
        None.

    Returns:
        None.
    """
    
    data_MPS = np.loadtxt('../Results/Matter_PS.txt')
    
    #k: wavenumber.
    
    k: np.array(float) = data_MPS[:, 0]
    
    #Matter_PS: matter power-spectrum.
    
    Matter_PS: np.array(float) = data_MPS[:, 1]
    
    #Read and assign the data from 'galaxy_survey_data.txt' in folder 'Data'.
    
    data_gal = np.loadtxt('../Data/galaxy_survey_data.txt', skiprows=1)
    
    k_gal: np.array(float) = data_gal[:, 0]
    
    P_gal: np.array(float) = data_gal[:, 1]
    
    ErrorP_gal: np.array(float) = data_gal[:, 2]
    
    #Read and assign the data from 'WMAP_ACT_data.txt' in folder 'Data'.
    
    data_ACT = np.loadtxt('../Data/WMAP_ACT_data.txt', skiprows=1)
    
    k_ACT: np.array(float) = data_ACT[:, 0]
    
    P_ACT: np.array(float) = data_ACT[:, 1]
    
    P_upper: np.array(float) = data_ACT[:, 2]
        
    #Read and assign the data from 'Ly_alpha.txt' in folder 'Data'.
    
    data_Ly = np.loadtxt('../Data/Ly_alpha.txt', skiprows=1)
    
    k_Ly: np.array(float) = data_Ly[:, 0]
    
    P_Ly: np.array(float) = data_Ly[:, 1]
    
    DeltaP_Ly: np.array(float) = data_Ly[:, 2]
        
    #Read the data from 'cosmology.txt' in folder 'Results'.
    
    data = np.loadtxt('../Results/cosmology.txt')
    
    #Get index of radiation and matter equality.
    
    index: int = aux.index_equality()[0]
    
    #Hp_eq: conformal Hubble factor when there is radiation and matter equality.
    
    Hp_eq: np.array(float) = data[index, 2]*10*sc.constants.parsec #100km/(Mpc*s)
    
    #k_eq: wavenumber when there is radiation and matter equality.
    
    k_eq: float = Hp_eq*1e5/(sc.constants.c*0.67) #h/Mpc
    
    #Make the plot.
        
    plt.figure()
    plt.plot(k, Matter_PS, label='Theory prediction')
    plt.errorbar(k_gal, P_gal, ErrorP_gal, label='SDSS Galaxies (DR7 LRG)',\
                 ls='none', marker='.', capsize=2)
    plt.errorbar(k_ACT, P_ACT, P_upper, label='CMB (WMAP+ACT)', ls='none',\
                 marker='.', capsize=2)
    plt.errorbar(k_Ly, P_Ly, DeltaP_Ly, label=r'Lyman $\alpha$ forest',\
                 ls='none', marker='.', capsize=2)
    plt.vlines(k_eq, plt.ylim()[0], plt.ylim()[1], label=r'$k_{eq}$',\
               ls='--', color='black')
    plt.xlabel(r'Wavenumber $k$ [$h$/Mpc]')
    plt.ylabel(r'$P(k)$ [(Mpc/$h$)$^3$]')
    plt.title('The total matter power-spectrum')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.grid()
    aux.show_or_save('../Plots/Milestone IV/Total matter PS.pdf')
    
def other_plots():
    """
    Generate additional plots related to cosmological analysis.

    This function reads data from a file containing information for plotting
    various cosmological quantities. It then generates plots to visualize the
    temperature multipoles and their squares divided by wavenumber.

    Parameters:
        None.
    Returns:
        None.
    """
    
    data = np.loadtxt('../Results/Matter_PS.txt')
    
    #k: wavenumber.
    
    k: np.array(float) = data[:, 0]
    
    #keta0: k times eta_0.
    
    keta0: np.array(float) = data[:, 2]
    
    #Theta_l: photon temperature multipoles.
    
    Theta_2: np.array(float) = data[:, 3]
    
    Theta_20: np.array(float) = data[:, 4]
    
    Theta_200: np.array(float) = data[:, 5]
    
    Theta_2000: np.array(float) = data[:, 6]
    
    plt.figure()
    plt.plot(keta0, Theta_2, label=r'$\Theta_2(k)$')
    plt.plot(keta0, Theta_20, label=r'$\Theta_{20}(k)$')
    plt.plot(keta0, Theta_200, label=r'$\Theta_{200}(k)$')
    plt.plot(keta0, Theta_2000, label=r'$\Theta_{2000}(k)$')
    plt.xlabel(r'$k\eta_0$')
    plt.ylabel(r'$\Theta_\ell$')
    plt.title(r'$\Theta_\ell$ for different $\ell$ values')
    plt.grid()
    plt.legend()
    aux.show_or_save('../Plots/Milestone IV/Theta_l.pdf')
    
    plt.figure()
    plt.plot(k, Theta_2**2/k, label=r'$|\Theta_2(k)|^2/k$')
    plt.plot(k, Theta_20**2/k, label=r'$|\Theta_{20}(k)|^2/k$')
    plt.plot(k, Theta_200**2/k, label=r'$|\Theta_{200}(k)|^2/k$')
    plt.plot(k, Theta_2000**2/k, label=r'$|\Theta_{2000}(k)|^2/k$')
    plt.xlabel(r'$k$ [m$^{-1}$]')
    plt.ylabel(r'$|\Theta_\ell(k)|^2/k$')
    plt.title(r'$\Theta_\ell^2/k$ for different $\ell$ values')
    plt.grid()
    plt.legend()
    aux.show_or_save('../Plots/Milestone IV/Theta_l_over_k.pdf')
    
def _planck_cmap():

    """
    Build the CMB colormap used in Planck's official temperature maps.

    If 'Data/Planck_Parchment_RGB.txt' is present (download it from
    https://github.com/zonca/paperplots/raw/master/data/Planck_Parchment_RGB.txt
    and place it in the Data/ folder) it is used directly: it is the exact
    color table used in the official Planck release. Otherwise, a close
    hand-built approximation of the same blue -> white/cream -> red
    "Parchment" diverging scheme is used instead, so the map always looks
    close to the Planck style even without the external file.

    Parameters:
        None.

    Returns:
        matplotlib.colors.Colormap
    """

    import os

    from matplotlib.colors import ListedColormap, LinearSegmentedColormap

    path = '../Data/Planck_Parchment_RGB.txt'

    if os.path.exists(path):

        cmap = ListedColormap(np.loadtxt(path)/255.)

    else:

        stops = [(0.00, '#0000A0'), (0.15, '#0000FF'), (0.30, '#00AAFF'),
                  (0.45, '#E6F5FF'), (0.50, '#FFFFF0'), (0.55, '#FFF5B0'),
                  (0.70, '#FFC000'), (0.85, '#FF4500'), (1.00, '#800000')]

        cmap = LinearSegmentedColormap.from_list('planck_like', stops)

    cmap.set_bad('gray')

    cmap.set_under('white')

    return cmap

def CMB_map():
    """
    Generate a map of the Cosmic Microwave Background (CMB).

    This function reads data from a file containing the CMB power spectrum
    information. It then generates a random realization of the CMB based on the
    power spectrum data, onverts spherical harmonic coefficients to a map, and
    plots the resulting CMB map.

    Returns:
        Parameters.
    Returns:
        None.
    """

    import healpy as hp
    
    nside = 2**10
    
    #Read the CMB power spectrum 'cells.txt' in folder 'Results'.
    
    data = np.loadtxt('../Results/cells.txt')
    
    #Extract the first column as l values
    
    ell = data[:, 0].astype(int)
    
    #Extract the second column as Cl values
    
    C_ell = data[:, 1]*(2*np.pi)/(ell*(ell+1))/((2.7255*1e6)**2)
    
    #'cells.txt' starts at ell=2 (there is no theoretical prediction for the
    #monopole or dipole), but healpy expects a 0-indexed array where entry i
    #is the power at multipole l=i, from l=0 up to l=lmax. Build that array
    #explicitly here (instead of passing C_ell directly, which would silently
    #get mis-read as starting from ell=0, shifting every multipole by 2).
    #Indices 0 and 1 (monopole and dipole) are left at zero: this is also
    #exactly what is done to the official Planck maps, since the monopole is
    #just the mean temperature and the observed dipole is dominated by our
    #own motion relative to the CMB rest frame, not by cosmological structure.
    
    ellmax = int(ell.max())
    
    Cl = np.zeros(ellmax+1)
    
    Cl[ell] = C_ell
    
    #Set the random seed to a specific value so the map is always the same
    
    np.random.seed(0)
    
    #Generate random spherical harmonic coefficients
    
    alm = hp.synalm(Cl, lmax=ellmax, new=True)
    
    #Convert the spherical harmonic coefficients to a map
    
    cmb_map = hp.alm2map(alm, nside)
    
    #'alm' (and therefore 'cmb_map') are in dimensionless Delta T/T units,
    #since Cl was divided by T_CMB^2 above (needed for synalm to behave).
    #Convert back to physical temperature units (muK, the standard unit for
    #CMB temperature maps) before plotting, otherwise the colorbar numbers
    #and the 'unit' label below would be meaningless.
    
    T_CMB_muK = 2.7255e6
    
    cmb_map_muK = cmb_map*T_CMB_muK
    
    hp.mollview(cmb_map_muK, title='The Cosmic Microwave Background',\
                cmap=_planck_cmap(), unit=r'$\mu K$', min=-300, max=300)
    
    aux.show_or_save('../Plots/Milestone IV/CMB map.pdf')
    
def milestone4():
    
    """
    Execute Milestone IV functionality.

    This function executes Milestone IV functionality, which involves plotting
    various cosmological quantities, including temperature multipoles, CMB
    power-spectrum, matter power-spectrum, and generating a map of the CMB.

    Parameters:
        None.

    Returns:
        None.
    """
    
    #Plot theta_l plots.
    
    other_plots()
    
    #Plot the CMB power-spectrum.
    
    CMB_PowerSpectrum()
    
    #Plot the different components, if they were computed (see
    #Constants.compute_PS_components in Source/Utils.h).
    
    if bool(aux.read_flags().get('compute_PS_components', 1)):
        
        CMB_PS_components()
    
    #Plot the matter power-spectrum.
    
    Matter_PowerSpectrum()

    #Plot the CMB map.
    
    CMB_map()
    
    #Show every plot from this milestone in a floating window, if requested
    #(CMB_SHOW_PLOTS=1); does nothing otherwise.
    
    aux.finalize_plots()
    
    #Let the user know (via Telegram) that the simulation has finished,
    #attaching the CMB TT power-spectrum as a quick visual summary. Wrapped
    #in try/except so that a network hiccup or missing internet access never
    #crashes the pipeline at the very last step.
    
    try:
        
        send_notification = bool(aux.read_flags().get('send_notification', 1))
        
        if not send_notification:
            
            pass
            
        elif notif is None:
            
            print('\nWarning: could not import notification.py (missing the "requests" package?). Notification skipped.')
            
        elif _last_cmb_ps_jpg is not None:
            
            notif.bot_texter('✅ Simulation finished.',
                              file_data=_last_cmb_ps_jpg,
                              file_name=_last_cmb_ps_filename)
            
        else:
            
            notif.bot_texter('✅ Simulation finished.')
            
    except Exception as e:
        
        print(f'\nWarning: could not send the Telegram notification ({e}).')

if __name__ == "__main__":
    
    #Run the main function of the script and execute the visualization of data.
    
    milestone4()
    