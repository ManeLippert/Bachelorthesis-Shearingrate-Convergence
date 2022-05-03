# --- BSC MareNostrum 4 ---
# Cores: 48 cores/node
# 1.88Gb/core (default nodes)
# 7.928Gb/core (high memory nodes)
# 
# Test cases ok, except:
#     chease_cf_modebox, last digit difference in one of the apar fluxes
#     global_em in parallel (code stalled, not sure why...)

REQUIREDMODULES= intel/2017.4  impi/2017.4  fftw/3.3.6  bsc/1.0  

MPI = usempi2
FFTLIB = FFT_FFTW3
#SMP = OPENMP

OPTFLAGS = on
DEBUG = off
PERF=perf

CC = icc
FC = mpiifort
LD = $(FC)

IMPLICIT=umfpack
#IO_LIB = HAVE_HDF5
GKW_LIBS=$GKW_HOME/libs

FFLAGS_DOUBLE = -r8
#FFLAGS_OPT = -no-prec-div # breaks geom related test cases
FFLAGS_OPT = -xCORE-AVX512 -mtune=skylake -O3 -ip -noalign
FFLAGS_OMP = -qopenmp
FFLAGS_INC = -I. -I/usr/include $(FFTW_INCL)
LDFLAGS    = -lfftw3 -lumfpack
FFLAGS_DEBUG = -fpe0 -O0 -g -traceback -check assume,pointers,stack,uninit,bounds

