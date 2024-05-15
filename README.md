# CMB power spectrum
## The Cosmic Microwave Background and the Large Scale Structure of our Universe
This repository will store the code used to compute the CMB power spectrum. I will update it as I progress with this project for the course AST5220 - Cosmology II of the University of Oslo (UiO). The structure of the repository and the code is based on Hans A. Whinter's (@HAWinther, https://cmb.wintherscoming.no/index.php). This repository contains all the needed files for runnig the project, however, only the code related to the second milestone has been completed. At the moment, the code can only solve the cosmological background, make a MCMC fit to the given supernova data and solve the recombinatioon history of the Universe.

Edit the Makefile adding the right paths and compile the code running [ make ]. The code runs from Main.cpp and then proceeds to go through the different milestones one by one untill we have the CMB power spectra in the end. If you get it compiled then run it as [ ./cmb ]. If you want to compile this on your computer you need to install the GSL library first.

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

- If this fails with "libgsl.so not found" then run the command: export LD\_LIBRARY\_PATH="$LD\_LIBRARY\_PATH:$HOME/local/lib"

and try to run ./cmb again and it should work. To avoid having to run this command every time you open a new terminal open the $HOME/.bashrc file and add this line to the end of the file and it will load everytime you open a new window.

## Plot the results

After running the code, you will need to run the Python codes, simply running Main.py. This will generate relevant plots and will print some relevant information.

## Read the reports

The folder named 'Reports' will store the reports of each milestone. At the moment, only the reports for Milestone I (I: Solving the cosmological background.pdf), Milestone II (II: Recombination history of the Universe.pdf) and Milestone III (III: Evolution of structure in the Universe.pdf) are available. Also a report summarizing the three articles is available (Computing the CMB power spectrum.pdf).
