# Curie Thin nodes  @ TGCC, 80640 cores (Intel Sandy Bridge 2.3GHz)


OPTFLAGS = on
DEBUG = off

FC = mpif90
FFLAGS_DEBUG = -g -traceback -C -ftrapuv
FFLAGS_OPT = -O3 -ipo -no-prec-div
FFLAGS_DOUBLE = -r8
FFLAGS_INC = -I$(FFTW3_INC_DIR)
FFLAGS_OMP = -openmp
LD = $(FC)
LDFLAGS = $(CCC_FFTW3_LDFLAGS) #-lfftw3f

REAL_PRECISION = real_precision_double
MPI = mpi2
FFTLIB = FFT_FFTW3
#IMPLICIT = umfpack
SMP = OPENMP

CC = cc
MPIRUNCMD = ccc_mprun

