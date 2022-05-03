## ---------- TOK cluster IPP Garching --------------
# 
MPI = usempi2
FFTLIB =FFT_FFTW3
REQUIREDMODULES = impi/4.0.0 pgi fftw-gcc
CC = cc
FC = mpipgf
FFLAGS_OPT = -O2
FFLAGS_DEBUG = -g -C -traceback
FFLAGS_DOUBLE = -r8
FFLAGS_INC = -I/usr/include -I$(FFTW_HOME)/include -I./
LDFLAGS = -L$(FFTW_HOME)/lib -lfftw3 -lfftw3f 

#comment the lines below if you don't need the global, nonspectral or implicit schemes
#LDFLAGS+=-L/afs/ipp/home/f/fjc/gkw/libs/UMFPACK/Lib/ -lumfpack_gkw64_NB -L/afs/ipp/home/f/fjc/gkw/libs/AMD/Lib/ -lamd
#IMPLICIT=umfpack

