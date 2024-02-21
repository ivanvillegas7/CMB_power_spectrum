# CMB power spectrum
This repository will store the code used to compute the CMB power spectrum. I will update it as I progress with this project for the course AST5220 - Cosmology II of the University of Oslo (UiO). The structure of the repository and the code is based on Hans A. Whinter's (HAWinther). This repository contains all the needed files for runnig the project, however, only the code related to the first milestone has been completed. At the moment, the project can only compute the cosmological background.

Compile the code running [ make ]. The code runs from Main.cpp and then proceeds to go through the different milestones one by one untill we have the CMB power spectra in the end. If you get it compiled then run it as [ ./cmb ]. If you want to compile this on your computer you need to install the GSL library first. See below for instructions if you haven't installed a library like this before.

On Linux or a Mac you can either use a package manager or install it directly as follows:

    Go the the home directory:

cd $HOME

    Make a local folder to hold libraries:

mkdir local

    Enter this directory:

cd local

    Download the code (if you don't have wget you need to get the file to this dir by other means):

wget ftp://ftp.gnu.org/gnu/gsl/gsl-2.6.tar.gz

    Untar the code:

tar -xvf gsl-2.6.tar.gz

    You should now have the gsl-2.6 folder. Enter it:

cd gsl-2.6

    Run the configure script:

./configure --prefix=$HOME/local

    Compile and install it:

make ; make install

    In the CMB code Makefile change the include and lib paths to point to the library:

INC = -I$(HOME)/local/include LIBS = -L$(HOME)/local/lib -lgsl -lgslcblas

    If this fails with "libgsl.so not found" then run the command:

export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/local/lib"

To avoid having to run this command every time you open a new terminal open the $HOME/.bashrc file and add this line to the end of the file and it will load everytime you open a new window.
