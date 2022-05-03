## -- Warwick AMD machine --
#REQUIREDMODULES = gnu/ompi64
#MPI = mpi2
FFTLIB = FFT_FFTW3
#SMP = OPENMP

CC = gcc
FC = mpif90
LD = $(FC)

FFLAGS_OPT = -O2 
FFLAGS_OMP = -fopenmp
FFLAGS_INC = -I/software/mathlib/gnu/include
LDFLAGS    = -L/usr/lib64 -lfftw3 -L/software/mathlib/gnu/lib -Wl,-rpath,/software/mathlib/gnu/lib 

