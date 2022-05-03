FC            = mpif90 
FFLAGS_WARN   = -W -Wall
FFLAGS_DEBUG  = -g -fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow
FFLAGS_OPT    = -O3
FFLAGS_DOUBLE = -fdefault-real-8 -fdefault-double-8
FFLAGS_INC    = -I.
FFLAGS_OMP    = -fopenmp
FFLAGS_OTHER  = -march=native

LD = $(FC)

GKW_LIBS=../../libs
LDFLAGS= -lfftw3 -lfftw3f $(GKW_LIBS)/UMFPACK/Lib/libumfpack_gkw64_NB.a $(GKW_LIBS)/AMD/Lib/libamd.a

REAL_PRECISION = real_precision_double
MPI            = usempi2
FFTLIB         = FFT_FFTW3
IMPLICIT       = umfpack
SMP            = OPENMP

CC = gcc

MPIRUNCMD = mpirun
