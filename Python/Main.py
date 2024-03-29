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
    
    print('Milestone I has been concluded.\n')
    
    Milestone2.milestone2()
    
    print('Milestone II has been concluded.\n')
    
    """
    Uncoment while completing milestones.
    """
    
    #Milestone3.milestone3()
    
    #print('Milestone III has been concluded.\n')
    
    #Milestone4.milestone4()
    
    #print('Milestone IV has been concluded.\n')
    
    #print('The project has been concluded.')
    
#Run the main function of the project and execute the visualization of data.

main()
