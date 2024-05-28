# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 12:11:52 2024

@author: Iván Villegas Pérez
"""

"""
This script contains the main function for the project. IT iports the files
related to each milestone and calls the main function of the milestone.
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
    
    Additionally, it can be uncommented to run subsequent milestones as they
    are completed.
    
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
    
    #Check if photon polarization has been included.
    
    polarization_txt: str
    
    polarization_txt = input('\nDid you include photon polarization? [Yes, no]: ')
    
    polarization: bool
    
    if polarization_txt.lower()=='no':
        
        polarization = False
        
    else:
        
        polarization = True
    
    Milestone3.milestone3(polarization)
    
    print('\nMilestone III has been concluded.\n')
    
    #Run the code related to computing the CMB and matter power-spectra.
    
    Milestone4.milestone4(polarization)
    
    print('Milestone IV has been concluded.\n')
    
    print('The project has been concluded.')
    
#Run the main function of the project and execute the visualization of data.

main()
