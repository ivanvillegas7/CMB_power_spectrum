# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 12:11:52 2024

@author: Iván Villegas Pérez
"""

"""
This script contains the main function for the project. It imports the files
related to each milestone and calls the main function of the milestone.

NOTE: since Source/Main.cpp now calls each Milestone script automatically
(via run_python()) right after the corresponding C++ module finishes, you
normally do NOT need to run this file by hand anymore. It is kept so the
whole pipeline can still be regenerated from Python alone (e.g. `python3
Main.py` from inside the Python/ folder) once Results/ has already been
populated by `./cmb`, without needing neutrinos/polarization to be typed in:
they are read automatically from 'Results/run_info.txt'.
"""

#Import all relevant files.

import Milestone1

import Milestone2

import Milestone3

import Milestone4
        
def main():
    
    """
    Main functionality for computing the CMB power spectrum.

    This function coordinates the execution of different milestones in the
    project. It runs the code related to solving the background cosmology of
    the Universe and shows the results related to supernova fitting.
    
    Parameters:
        None.

    Returns:
        None.
    """
    
    #Run the code related to solving the background cosmology of the Universe
    #and show the results related to the supernove fitting.
    
    Milestone1.milestone1()
    
    #Print a message saying that the execution of the code is over.
    
    print('\nMilestone I has been concluded.')
    
    #Run the code related to solving the recombination history of the Universe.
    
    Milestone2.milestone2()
    
    print('\nMilestone II has been concluded.\n')
    
    #Run the code related to solving the perturbations of the Universe.
    #Neutrinos/polarization flags are read automatically from
    #'Results/run_info.txt' inside Milestone3, no user input needed.
    
    Milestone3.milestone3()
    
    print('\nMilestone III has been concluded.\n')
    
    #Run the code related to computing the CMB and matter power-spectra.
    #Same as above: polarization is read automatically from file.
    
    Milestone4.milestone4()
    
    print('Milestone IV has been concluded.\n')
    
    print('The project has been concluded.')
    
if __name__ == "__main__":
    
    #Run the main function of the project and execute the visualization of data.
    
    main()
