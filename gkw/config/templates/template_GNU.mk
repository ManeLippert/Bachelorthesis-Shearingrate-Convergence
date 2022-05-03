## -- GNU compiler template --
#REQUIREDMODULES= GNU compiler, MPI

MPI = usempi2
FFTLIB = FFT_FFTW3
#SMP = OPENMP
#DEBUG = on

FC = mpif90
FFLAGS_INC = -I. -I/usr/include -I${FFTW_HOME}/include
FFLAGS_OTHER = -std=f95 #-W -Wall -Wunderflow
FFLAGS_OMP = -fopenmp
FFLAGS_OPT = -O3 -ftracer -fomit-frame-pointer -pipe -fweb
FFLAGS_DEBUG = -g -fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow
FFLAGS_DOUBLE = -fdefault-real-8 -fdefault-double-8
LDFLAGS = -L${FFTW_HOME}/lib -lfftw3 -lfftw3f

# Make code timings, no reason not to use
PERF=perf

# Need the lines below if you want to use nonspectral or implicit scheme
# GKW_LIBS = /path/to/compiled/gkw/libs   # compile them with the same compiler 
# IMPLICIT=umfpack
# LDFLAGS+=-L${GKW_LIBS}/UMFPACK/Lib/ -lumfpack_gkw64_NB -L${GKW_LIBS}/AMD/Lib/ -lamd
