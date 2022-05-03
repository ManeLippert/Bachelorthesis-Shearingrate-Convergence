## -- Warwick cluster --

MODULES = intel/ompi64
MPI = mpi2
FFTLIB = FFT_FFTW3
#SMP = OPENMP

CC = cc
FC = mpif90
LD = $(FC)

FFLAGS_DOUBLE = -r8
FFLAGS_OPT = -O3 -no-prec-div
FFLAGS_OTHER = -vec-report0
FFLAGS_OMP = -openmp
FFLAGS_INC = -I/software/mathlib/intel/include
LDFLAGS    = -L/usr/lib64 -lpthread -liomp5 -lfftw3 -lfftw3f
FFLAGS_DEBUG = -g -C -traceback
