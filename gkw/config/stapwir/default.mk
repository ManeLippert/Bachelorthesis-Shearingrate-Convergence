OPTFLAGS = on
DEBUG = off
WARN = off

FC = mpif90
FFLAGS_DEBUG = -g -fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow
FFLAGS_OPT = -O3
FFLAGS_DOUBLE = -fdefault-real-8 -fdefault-double-8
FFLAGS_INC = -I. -I/usr/include
FFLAGS_OMP = -fopenmp
FFLAGS_OTHER = -march=native
FFLAGS_WARN  = -W -Wall
LD = $(FC)
FFLAGS_LD =
CC     = gcc
CFLAGS = -Ofast

REAL_PRECISION = real_precision_double
MPI            = usempi2
FFTLIB         = FFT_FFTW3
IMPLICIT       = umfpack
#SLEPC         = HAVE_SLEPC
#IO_LIB        = HAVE_HDF5
SMP            = OPENMP
#PERF          = perf
MPIRUNCMD      = mpirun
#INPUT_CHECK   = input_check
#UMFPACK=external
LD_UMFPACK=-lumfpack_gkw64_NB  -lcholmod  -lamd -lsuitesparseconfig -lcolamd -lblas
LD_BLAS=-lblas
LDFLAGS += -lfftw3 -lfftw3f
