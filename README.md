# CMB power spectrum
This repository will store the code used to compute the CMB power spectrum. I will update it as I progress with this project for the course AST5220 - Cosmology II of the University of Oslo (UiO). The structure of the repository and the code is based on Hans A. Whinter's (@HAWinther, https://cmb.wintherscoming.no/index.php). This repository contains all the needed files for runnig the project, however, only the code related to the first milestone has been completed. At the moment, the code can only solve the cosmological background and make a MCMC fit to the given supernova data.

Edit the Makefile adding the right paths and compile the code running [ make ]. The code runs from Main.cpp and then proceeds to go through the different milestones one by one untill we have the CMB power spectra in the end. If you get it compiled then run it as [ ./cmb ]. If you want to compile this on your computer you need to install the GSL library first.

After running the code, you will need to run the Python codes, simply running Main.py. This will generate relevant plots and will print some relevant information.
