OPTFLAGS = on
DEBUG    = off

#SMP      = off
MPI = mpi2
FFTLIB   = FFT_FFTW3
#HDF5=HDF5
#IMPLICIT=umfpack

CC = cc

ifeq ($(FFTLIB),FFT_FFTW3)
  LDFLAGS       = -lfftw3 -lfftw3f
endif

COMPILER = gnu
