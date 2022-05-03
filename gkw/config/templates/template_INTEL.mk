## -- Intel compiler template
#REQUIREDMODULES= Intel compiler, MPI, FFTW

MPI = usempi2
FFTLIB = FFT_FFTW3

#SMP=OPENMP
#DEBUG=on

CC = cc
FC = mpif90  # or mpiifort on many machines
LD = $(FC)
FFLAGS_OPT = -O3 -no-prec-div
FFLAGS_DEBUG = -traceback -g -C -ftrapuv
FFLAGS_DOUBLE = -r8
FFLAGS_OMP = -openmp
FFLAGS_INC = -I/usr/include -I./ #-L${FFTW_HOME}/include
FFLAGS_OTHER = -vec-report0 -std95
LDFLAGS = -L/usr/lib -lfftw3 -lfftw3f #-L${FFTW_HOME}/lib

# Make code timings, no reason not to use
PERF=perf

# Need the lines below if you want to use nonspectral or implicit scheme
# GKW_LIBS = /path/to/compiled/gkw/libs   # compile them with the same compiler 
# IMPLICIT=umfpack
# LDFLAGS+=-L${GKW_LIBS}/UMFPACK/Lib/ -lumfpack_gkw64_NB -L${GKW_LIBS}/AMD/Lib/ -lamd
