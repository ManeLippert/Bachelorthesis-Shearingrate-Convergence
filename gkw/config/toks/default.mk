## ---------- TOK cluster IPP Garching --------------
# 
MPI = usempi2
FFTLIB =FFT_FFTW3
REQUIREDMODULES = impi/4.0.0  intel/11.1 fftw/3.2.2
CC = cc
FC = mpiifort
#FFLAGS_OPT = -O3 -no-prec-div # These optimisations break the miller geometry test case with intel compiler 12.1
FFLAGS_OPT = -O2
FFLAGS_DEBUG = -g -C -traceback
FFLAGS_DOUBLE = -r8
FFLAGS_OMP = -openmp
FFLAGS_INC = -I/usr/include -I$(FFTW_HOME)/include -I./
FFLAGS_OTHER = -vec-report0 -std95
LDFLAGS = -L$(FFTW_HOME)/lib -lfftw3 -lfftw3f 
#SMP=OPENMP

#comment the lines below if you don't need the global, nonspectral or implicit schemes
LDFLAGS+=-L/afs/ipp/home/f/fjc/gkw/libs/UMFPACK/Lib/ -lumfpack_gkw64_NB -L/afs/ipp/home/f/fjc/gkw/libs/AMD/Lib/ -lamd
IMPLICIT=umfpack

