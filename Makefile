# Hans A. Winther (2020) (hans.a.winther@gmail.com)

SHELL:= /bin/bash

# Set compiler (use =c++17 if you have this availiable)
CC = g++ -std=c++17 

# Paths to GSL library (MIGHT NEED TO BE CHANGED)
INC  = -I/usr/local/include
LIBS = -L/usr/local/lib -lgsl -lgslcblas

#=======================================================
# Options
#=======================================================

# Show warnings if atempting to evaluate a spline out of bounds
COMMON_OPTIONS = -D_SPLINE_WARNINGS_ON

# Show info about the solution as we integrate
# COMMON_OPTIONS += -D_FIDUCIAL_VERBOSE_ODE_SOLVER_TRUE

# Add OpenMP parallelization. The perturbation ODE solve (loop over k) and
# the line-of-sight integration (also looped over k, now done up to 5 times
# per run for the TT/SW/ISW/DOPPLER/POL/E components) are the most expensive
# steps in the whole pipeline, and are already written as
# "#pragma omp parallel for" loops in Perturbations.cpp / PowerSpectrum.cpp.
# Without -fopenmp the compiler just ignores those pragmas and runs single
# threaded. This does not change a single number in the output, only how
# many CPU cores are used to get there.
COMMON_OPTIONS += -D_USEOPEMP
CC += -fopenmp

#=======================================================
# Release build (default): use this to actually run the pipeline and
# generate the results/plots for the report. No bounds-checking overhead.
#=======================================================
RELEASE_OPTIONS = $(COMMON_OPTIONS)
C_RELEASE = -O3 -march=native -DNDEBUG $(RELEASE_OPTIONS)

#=======================================================
# Debug build: use this (make debug) while developing/chasing bugs. Adds
# STL bounds-checking (catches out-of-range vector access, invalid
# iterators, etc., like the indexing bugs fixed in Perturbations.cpp) at
# the cost of being noticeably slower. Keeps debug symbols (-g).
#=======================================================
DEBUG_OPTIONS = $(COMMON_OPTIONS) -D_GLIBCXX_DEBUG
C_DEBUG = -O0 -g $(DEBUG_OPTIONS)

#=======================================================

VPATH=Source/
TARGETS := cmb
all: $(TARGETS)

# OBJECT FILES (release, in ./obj_release)
OBJS         = obj_release/Main.o obj_release/Utils.o obj_release/BackgroundCosmology.o obj_release/RecombinationHistory.o obj_release/Perturbations.o obj_release/PowerSpectrum.o obj_release/Spline.o obj_release/ODESolver.o
# OBJECT FILES (debug, in ./obj_debug)
OBJS_DEBUG   = obj_debug/Main.o obj_debug/Utils.o obj_debug/BackgroundCosmology.o obj_debug/RecombinationHistory.o obj_debug/Perturbations.o obj_debug/PowerSpectrum.o obj_debug/Spline.o obj_debug/ODESolver.o

# DEPENDENCIES
Main.o                  : BackgroundCosmology.h RecombinationHistory.h Perturbations.h PowerSpectrum.h
Spline.o                : Spline.h
ODESolver.o             : ODESolver.h
Utils.o                 : Utils.h Spline.h ODESolver.h
BackgroundCosmology.o   : BackgroundCosmology.h Utils.h Spline.h ODESolver.h
RecombinationHistory.o  : RecombinationHistory.h BackgroundCosmology.h
Perturbations.o         : Perturbations.h BackgroundCosmology.h RecombinationHistory.h
PowerSpectrum.o         : PowerSpectrum.h BackgroundCosmology.h RecombinationHistory.h Perturbations.h
Examples.o              : Utils.h Spline.h ODESolver.h

examples: Examples.o Utils.o Spline.o ODESolver.o
	${CC} -o $@ $^ $(C_RELEASE) $(INC) $(LIBS)

# Fast build used for actually generating results (make / make cmb / make all)
cmb: $(OBJS)
	${CC} -o $@ $^ $(C_RELEASE) $(INC) $(LIBS)

obj_release/%.o: %.cpp | obj_release
	${CC} -c -o $@ $< $(C_RELEASE) $(INC)

obj_release:
	mkdir -p obj_release

# Bounds-checked build used for debugging (make debug)
debug: $(OBJS_DEBUG)
	${CC} -o cmb_debug $^ $(C_DEBUG) $(INC) $(LIBS)

obj_debug/%.o: %.cpp | obj_debug
	${CC} -c -o $@ $< $(C_DEBUG) $(INC)

obj_debug:
	mkdir -p obj_debug

clean:
	rm -rf $(TARGETS) cmb_debug obj_release obj_debug *.o

