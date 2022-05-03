## -- CCFE CULHAM cluster, intel compiler  --
# This non longer works becuase GKW now requires mpi2
# This machine has not had any software updates since 2009 !

MPI = mpi
FFTLIB = FFT_FFTW3
REQUIREDMODULES = mpich/mx/64bit/INTEL mx-driver/64bit intel/compiler91_x86_64
CC = cc
FC = mpif90
FFLAGS_OPT = -O2
FFLAGS_DOUBLE = -r8
FFLAGS_INC = -I/usr/local/fusion/64/include
LDFLAGS = -L/usr/local/fusion/64/lib -lfftw3 -lfftw3f
