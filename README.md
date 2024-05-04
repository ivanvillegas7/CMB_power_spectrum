# CMB power spectrum
This repository will store the code used to compute the CMB power spectrum. I will update it as I progress with this project for the course AST5220 - Cosmology II of the University of Oslo (UiO). The structure of the repository and the code is based on Hans A. Whinter's (@HAWinther, https://cmb.wintherscoming.no/index.php). This repository contains all the needed files for runnig the project, however, only the code related to the second milestone has been completed. At the moment, the code can only solve the cosmological background, make a MCMC fit to the given supernova data and solve the recombinatioon history of the Universe.

Edit the Makefile adding the right paths and compile the code running [ make ]. The code runs from Main.cpp and then proceeds to go through the different milestones one by one untill we have the CMB power spectra in the end. If you get it compiled then run it as [ ./cmb ]. If you want to compile this on your computer you need to install the GSL library first.

After running the code, you will need to run the Python codes, simply running Main.py. This will generate relevant plots and will print some relevant information.

The folder named 'Reports' will store the reports of each milestone. At the moment, only the reports for Milestone I (I: Solving the cosmological background.pdf), Milestone II (II: Recombination history of the Universe.pdf) and Milestone III (III: Evolution of structure in the Universe.pdf) are available. Also a report summarizing the three articles is available (Computing the CMB power spectrum.pdf).
