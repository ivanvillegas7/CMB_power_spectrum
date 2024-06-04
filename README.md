# Computing the CMB power spectrum

## The Cosmic Microwave Background and the Large Scale Structure of our Universe

This repository will store the code used to compute the CMB power spectrum, a project for the course AST5220 - Cosmology II of the University of Oslo (UiO). The structure of the repository and the code is based on Hans A. Whinter's (@HAWinther, https://cmb.wintherscoming.no/index.php). This repository contains all the needed files for runnig the project, solving the cosmological background, making a MCMC fit to the given supernova data, solving the recombination history and solving the evolution of the Universe, as well as computing the CMB power-spectrum, the matter power-spectrum and generating a CMB map.

At the moment, the code is only available at the Master's degree version, i.e., it does not include neutrinos, Helium, reionization nor photon polarization. I will try to upgrade to the PhD version as soon as I have time to retake this project.

## Compiling and running

Edit the Makefile adding the right paths and compile the code running [ make ]. The code runs from Main.cpp and then proceeds to go through the different milestones one by one untill we have the CMB power spectra in the end. If you get it compiled then run it as [ ./cmb ]. If you want to compile this on your computer you need to install the [GSL library](ftp://ftp.gnu.org/gnu/gsl/) first.

## How to install GSL

See [this](https://solarianprogrammer.com/) for how to install it on a Windows machine. On Linux or a Mac you can either use a package manager or install it directly as follows:

- Go the the home directory: cd $HOME

- Make a local folder to hold libraries: mkdir local

- Enter this directory: cd local

- Download the code (if you don't have wget you need to get the file to this dir by other means): wget ftp://ftp.gnu.org/gnu/gsl/gsl-2.6.tar.gz

- Untar the code: tar -xvf gsl-2.6.tar.gz

- You should now have the gsl-2.6 folder. Enter it: cd gsl-2.6

- Run the configure script: ./configure --prefix=$HOME/local

- Compile and install it: make ; make install

- In the CMB code Makefile change the include and lib paths to point to the library:

  -INC  = -I$(HOME)/local/include
  
  -LIBS = -L$(HOME)/local/lib -lgsl -lgslcblas

- If this fails with "libgsl.so not found" then run the command:
  export LD\_LIBRARY\_PATH="$LD\_LIBRARY\_PATH:$HOME/local/lib"

and try to run ./cmb again and it should work. To avoid having to run this command every time you open a new terminal open the $HOME/.bashrc file and add this line to the end of the file and it will load everytime you open a new window.

## Plot the results

After running the code, you will need to run the Python codes, simply running Main.py. This will generate relevant plots and will print some relevant information. If you want to plot the results, you will need to install 'healpy'. You can do it simply running "conda install healpy" in your favourite conda distribution.

## Further documentation

The folder named 'Reports' will store the reports of each milestone. All the reports, for Milestone I (I: Solving the cosmological background.pdf), Milestone II (II: Recombination history of the Universe.pdf), Milestone III (III: Evolution of structure in the Universe.pdf) and Milestone IV (IV: The CMB and matter power-spectra.pdf), are available. Also an article summarizing the four reports is available (Computing the CMB power spectrum.pdf).
