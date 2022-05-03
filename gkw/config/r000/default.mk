 ## --- CINECA Marconi ---
#Marconi is the new Tier-0 system, co-designed by Cineca and based on
#the Lenovo NeXtScale platform, that substitutes the former IBM BG/Q
#system (FERMI) . MARCONI is based on the next-generation of the Intel®
#Xeon Phi™ product family alongside with Intel® Xeon® processor
#E5-2600 v4 product family.
# Model: Lenovo NeXtScale
# Architecture: Intel OmniPath Cluster
# Nodes: 1.512
# Processors: 2 x 18-cores Intel Xeon E5-2697 v4 (Broadwell) at 2.30 GHz
# Cores: 36 cores/node, 54.432 cores in total
# RAM: 128 GB/node, 3.5 GB/core
# Internal Network: Intel OmniPath

OPTFLAGS = on 
DEBUG = off 
PERF=perf

#these modules are recommended, not enforced. Load these in your ~/.profile or thelike.
REQUIREDMODULES= env-skl intel/pe-xe-2018--binary profile/advanced autoload intel intelmpi fftw petsc/3.8.3_complex--intelmpi--2018--binary slep
c/3.8.2_complex--intelmpi--2018--binary zlib/1.2.11--intel--pe-xe-2017--binary szip hdf5/1.8.17--intel--pe-xe-2017--binary

# These are for the eigenvalue solver
#SLEPC = HAVE_SLEPC
#FFLAGS_INC_SLEPC= -I$(PETSC_INCLUDE) -I$(SLEPC_INCLUDE)
#LDFLAGS_SLEPC=-Wl,-rpath,${SLEPC_LIB} -L${SLEPC_LIB} -lslepc -Wl,-rpath,${PETSC_LIB} -L${PETSC_LIB} -lpetsc #-lflapack -lfblas

# These are for the HDF5 format:
# IO_LIB = HAVE_HDF5
# The compiler flags for HDF5 are taken from the wrapper provided by HDF5
# FFLAGS_INC_H5= $(shell h5fc -show | sed 's/^[^ ]*//')
# LDFLAGS_H5=$(shell h5fc -show | sed 's/^[^ ]*//')
FC=mpiifort


# These are for OpenMP
SMP = OPENMP
FFLAGS_OMP = -fopenmp

# Flags and paths to use FFTW3 as the fft library:
FFTLIB = FFT_FFTW3
FFLAGS_INC_FFT= -I$(FFTW_INCLUDE)
LDFLAGS_FFT=-L$(FFTW_LIB) -lfftw3

# These are for umfpack
IMPLICIT=umfpack
LDFLAGS_UMFPACK=-lumfpack

# These are for the MPI library
MPI = usempi2


FFLAGS_DEBUG = -fpe0 -O0 -g -traceback -check assume,pointers,stack,uninit,bounds
FFLAGS_OPT = -xCORE-AVX512 -mtune=skylake -O3 -ip 
#FFLAGS_OPT = -02 -ip #-no-prec-div #-xAVX #problems in umfpack, and no faster
FFLAGS_DOUBLE = -r8

FFLAGS_INC = -I. $(FFLAGS_INC_SLEPC) $(FFLAGS_INC_FFT) $(FFLAGS_INC_H5)
LDFLAGS = $(LDFLAGS_SLEPC) $(LDFLAGS_FFT) $(LDFLAGS_UMFPACK) $(LDFLAGS_H5)
