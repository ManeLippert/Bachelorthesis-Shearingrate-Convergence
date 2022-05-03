## -- Warwick desktop FJC --
REQUIREDMODULES = intel/11.1 intel/mpich2/1.2.1p1
MPI = mpi2
FFTLIB = FFT_FFTW3
SMP = OPENMP

CC = cc
FC = mpif90
LD = $(FC)
FFLAGS_OPT = -O3 -no-prec-div
FFLAGS_DEBUG = -traceback -g -C -ftrapuv
FFLAGS_DOUBLE = -r8
FFLAGS_OMP = -openmp
FFLAGS_INC =
FFLAGS_OTHER = -vec-report0 
LDFLAGS = -L/warwick/intel/Compiler/11.1/073/lib/intel64 -lfftw3 -lfftw3f
