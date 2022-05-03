# -- intel atom netbook GNU/Linux --

## default compiler
COMPILER = gnu

OPTFLAGS = on
DEBUG = on

REAL_PRECISION = real_precision_double
MPI = mpi2
FFTLIB = FFT_FFTW3
SMP = OPENMP

MPIRUNCMD = mpirun
CC = gcc
EXEC_PREFIX = gkw
EXEC_SUFFIX = .x
