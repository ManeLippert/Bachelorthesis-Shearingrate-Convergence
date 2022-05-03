# Bayreuth new cluster - gfortran (mpif90) and SLEPC, no FFTW
 
## Without MPI, only single processor runs can be performed 
## Uncomment to use MPI:
MPI = mpi2

## Without FFTW (default), only linear runs can be performed (off by default) 
## Uncomment to use FFTW:
#FFTLIB = FFT_FFTW3

## Directory where the fftw library is located. 
# This one is for GNU compiler ?
#FFT_LIB=/home/vivaorg/vivaorg_10150057/fftw-3.3.3

# This one is for intel compiler ?
#FFT_LIB=/home/40/bt121140/fftw

# Defining directory of slepc and petsc
PETSC_DIR=/home/vivaorg/vivaorg_10150057/local_mpif90
SLEPC_DIR=/home/vivaorg/vivaorg_10150057/local_mpif90

## Turn compiler debugging on (off by default)
DEBUG = off

## Turn compiler optismations off (on by default)
OPTFLAGS = on

## Use the gnu compiler settings in subfolder as the default
#COMPILER=gnu

## Linking flags ( e.g. -L/path/to/include -lfftw3)
LDFLAGS = -L${PETSC_DIR}/lib -lslepc -lpetsc -lflapack -lfblas \
          -L../../libs/UMFPACK/Lib/ -lumfpack_gkw64_NB -L../../libs/AMD/Lib/ -lamd

## To compile in single precision, use real_precision_default
#REAL_PRECISION = real_precision_default

# For machines which require enviroment modules to be loaded, 
# Optionally list the modules that should be loaded for use with this config
# Make will check and report if correct modules are not loaded.
# The user is required to load the correct modules before compilation
REQUIREDMODULES = intel-cluster-studio-2013 \
                   mpi/intel/4.1.0.030

#SMP = OPENMP
IMPLICIT = umfpack

SLEPC = HAVE_SLEPC

## Override name of linker (default is to use fortran compiler)
#LD

##Name of c compiler (for svnrev version tracker)
#CC = cc

## Name for GKW if the version number cannot be found with svnrev
#DEFAULT_GKW_VERSION = SOME_VERSION

#PERF=perf

CC = gcc
FC = /cluster/intel/impi/4.1.0.030/bin64/mpif90
LD = $(FC)

GKW_LIBS=$HOME/gkw/trunk/libs

FFLAGS_DOUBLE = -fdefault-real-8 -fdefault-double-8
FFLAGS_OPT = -O3
FFLAGS_DEBUG = -g
FFLAGS_OTHER = 
FFLAGS_OMP = -fopenmp

## Compiler include flags ( e.g. -I/path/to/include)
FFLAGS_INC = -I. -I/cluster/intel/impi/4.1.0.030/include64 -I$(PETSC_DIR)/include
