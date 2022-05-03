## -- Bayreuth cluster --

REQUIREDMODULES =  mpi/openmpi/1.3.2/intel
MPI = mpi2
FFTLIB = FFT_FFTW3
IMPLICIT = umfpack
#SMP = OPENMP

CC = cc
FC = mpif90
LD = $(FC)

FFLAGS_DOUBLE = -r8
FFLAGS_OPT = -O3 -no-prec-div
FFLAGS_OTHER = -vec-report0 #-Warn all
FFLAGS_OMP = -openmp
FFLAGS_INC = -I/software/mathlib/intel/include
LDFLAGS    = -L/usr/lib64 -lpthread -lfftw3 -L../../libs/UMFPACK/Lib -lumfpack_gkw64_NB -L../../libs/AMD/Lib -lamd
FFLAGS_DEBUG = -g -C -traceback
