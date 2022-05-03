## -- GNU compiler on tok cluster --
#REQUIREDMODULES= GNU compiler, MPI

MPI = mpi2
FFTLIB = FFT_FFTW3
#SMP = OPENMP
#DEBUG = on

FFTW_HOME=/afs/@cell/common/soft/fftw/fftw-3.3.3/@sys/gcc-4.7/impi-4.1

FC = mpif90
FFLAGS_INC = -I. -I/usr/include -I${FFTW_HOME}/include
#FFLAGS_OTHER = -std=f95 -pedantic-errors # (works only with IMPLICIT=noimp)
#FFLAGS_OTHER = -W -Wall -Wunderflow 
FFLAGS_OMP = -fopenmp
FFLAGS_OPT = -O3 -ftracer -fomit-frame-pointer -pipe -fweb
FFLAGS_DEBUG = -g -fbounds-check -fbacktrace -fcheck=mem -fcheck=do #-ffpe-trap=zero,overflow,underflow
FFLAGS_DOUBLE = -fdefault-real-8 -fdefault-double-8
LDFLAGS = -L${FFTW_HOME}/lib -lfftw3 -lfftw3f

# Make code timings, no reason not to use
PERF=perf

# Need the lines below if you want to use nonspectral or implicit scheme

GKW_LIBS=/afs/ipp/home/f/fjc/gkw_clean/libs
# these compiled with the GNU compiler 

IMPLICIT=umfpack
LDFLAGS+=-L${GKW_LIBS}/UMFPACK/Lib/ -lumfpack_gkw64_NB -L${GKW_LIBS}/AMD/Lib/ -lamd
