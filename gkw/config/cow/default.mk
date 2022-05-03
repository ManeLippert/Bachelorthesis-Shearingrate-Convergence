## -- Warwick desktop cow, intel compiler
REQUIREDMODULES = intel/11.1 intel/mpich2/1.2.1p1
MPI = mpi2
FFTLIB = FFT_FFTW3

CC = cc
FC = mpif90
LD = $(FC)
FFLAGS_OPT = -O2 -static
FFLAGS_DEBUG = -traceback -g -C -ftrapuv
FFLAGS_DOUBLE = -r8
FFLAGS_INC =
FFLAGS_OTHER = -vec-report0
LDFLAGS = -lfftw3 -lfftw3f
