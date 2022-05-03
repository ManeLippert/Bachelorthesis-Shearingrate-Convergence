## -- Cullix desktop --
REQUIREDMODULES= openmpi-i386
# To get the FFTW include files you need to do the following:
# sudo yum install fftw-devel.i686

MPI = mpi2
FFTLIB = FFT_FFTW3
SMP = OPENMP
#DEBUG = on

FC = mpif90
FFLAGS_INC = -I. -I/usr/include 
FFLAGS_OTHER = -std=f95
FFLAGS_WARN = -fmax-errors=2 #-Wunused #-W -Wall -Wunderflow
FFLAGS_OMP = -fopenmp
FFLAGS_OPT = -O2 #-O3 -ftracer -fomit-frame-pointer -pipe -fweb
FFLAGS_DEBUG = -g -fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow
FFLAGS_DOUBLE = -fdefault-real-8 -fdefault-double-8
LDFLAGS = -L/usr/lib -lfftw3 -lfftw3f

# Make code timings, no reason not to use
PERF=perf

# Need the lines below if you want to use nonspectral or implicit scheme
GKW_LIBS=${GKW_HOME}/libs
IMPLICIT=umfpack32
#LDFLAGS+=-L${GKW_LIBS}/UMFPACK/Lib/ -lumfpack_gkw32_NB -L${GKW_LIBS}/AMD/Lib/ -lamd
