## -- Archer HPC UK resource: defaults --

# Do this before compiling
#module swap PrgEnv-cray PrgEnv-intel
#module load fftw

REQUIREDMODULES+=fftw
OPTFLAGS=on
SMP=OPENMP 

MPI = mpi2
FFTLIB = FFT_FFTW3
CC=gcc
PERF=perf

IMPLICIT=umfpack
GKW_LIBS=/home/e281/e281/phrhac/gkw/libs
LDFLAGS += -L${GKW_LIBS}/UMFPACK/Lib/ -lumfpack_gkw64_NB -L${GKW_LIBS}/AMD/Lib/ -lamd

#COMPILER = cray# cray fastest, and matches default loaded modules (but has io record length problems)
COMPILER=intel
