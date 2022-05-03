REQUIREDMODULES=gcc/5.0

FC = h5pfc
FFLAGS_WARN   = -W -Wall
FFLAGS_DEBUG  = -g -fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow
FFLAGS_OPT    = -Ofast
FFLAGS_DOUBLE = -fdefault-real-8 -fdefault-double-8
FFLAGS_INC    = -I. -I/opt/gcc-5.0/include
FFLAGS_OMP    = -fopenmp
FFLAGS_OTHER  = -march=bdver2

LD = $(FC)

GKW_LIBS=../../libs
LDFLAGS= -L/opt/gcc-5.0/lib -lfftw3 -lfftw3f $(GKW_LIBS)/UMFPACK/Lib/libumfpack_gkw64_NB.a $(GKW_LIBS)/AMD/Lib/libamd.a

REAL_PRECISION = real_precision_double
MPI            = usempi2
FFTLIB         = FFT_FFTW3
IMPLICIT       = umfpack
#SLEPC         = HAVE_SLEPC
IO_LIB         = HAVE_HDF5
SMP            = OPENMP

CC = gcc

MPIRUNCMD = mpirun
