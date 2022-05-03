## -- Yann's machine @ Marseille --
# Note: need to use mpich, test cases in parallel stall with openmpi
# And configure mpich with --with-device=ch3:sock to avoid stalling when cores are oversubscribed

MPI = usempi2
FFTLIB = FFT_FFTW3
OPTFLAGS = on
DEBUG = off
PERF=perf

CC = gcc
FC = mpifort
LD = $(FC)

IMPLICIT=umfpack
#IO_LIB = HAVE_HDF5
GKW_LIBS=/home/yann/codes/gkw/libs

FFLAGS_DOUBLE = -fdefault-real-8 -fdefault-double-8
FFLAGS_OPT = -O3 -ftracer -fomit-frame-pointer -pipe -fweb
FFLAGS_OTHER = #-W -Wall -Wunderflow
FFLAGS_OMP = -fopenmp
FFLAGS_INC = -I. -I/usr/include -I/usr/local/include
LDFLAGS    = -lfftw3 -lfftw3f -L$(GKW_LIBS)/UMFPACK/Lib -lumfpack_gkw64_NB -L$(GKW_LIBS)/AMD/Lib -lamd
FFLAGS_DEBUG = -g -Og -Wall -fimplicit-none -fcheck=all -fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow


##############################################################################
### GKW build configuration template                                       ###
##############################################################################

## Remove this line to avoid warnings AFTER the rest has been set up
## correctly and the code can be built.
#CONFIG_WARNING=yes

#==============================================================================
## Common configuration options, uncomment and edit those that are required.
#==============================================================================

## MPI and FFTW options require external libraries

## Without MPI, only single processor runs can be performed 
## Uncomment to use MPI:
#MPI = mpi2

## Without FFTW (default), only linear runs can be performed (off by default) 
## Uncomment to use FFTW:
#FFTLIB = FFT_FFTW3

## Turn compiler debugging on (off by default)
#DEBUG = on

## Turn compiler optismations off (on by default)
#OPTFLAGS = off

## Name of fortran compiler (default is gfortran)
#FC = mpif90

## Compiler debugging flags used if DEBUG=on (e.g. -g -C -traceback)
#FFLAGS_DEBUG = 

## Compiler optimisation flasgs used if OPTFLAGS=on (e.g -03)
#FFLAGS_OPT = -03

## Compiler flag used when compiling in double precision (default is for gfortran)
#FFLAGS_DOUBLE= -r8

## If not using a wrapper (such as h5pfc or mpif90), 
## appropriate FFLAGS_INC and LDFLAGS for the libraries may need to be set.

## Compiler include flags ( e.g. -I/path/to/include)
#FFLAGS_INC =

## Linking flags ( e.g. -L/path/to/include -lfftw3)
#LDFLAGS = 

## To compile in single precision, use real_precision_default
#REAL_PRECISION = real_precision_default

## For machines which require enviroment modules to be loaded, 
## Optionally list the modules that should be loaded for use with this config
## Make will check and report if correct modules are not loaded.
## The user is required to load the correct modules before compilation
#REQUIREDMODULES = 

##===========================================================================
## You are less likely to need to change anything below here:
##===========================================================================
## More preprocessor options unlikely to be needed in current standard usage
#SMP = OPENMP
#IMPLICIT = umfpack
#HDF5 = HDF5
#PERF = perf  # performance timings
#PREPROC_FLAG
#PREPROC_SEP
#PREPROC_PREFIX

## Flags required for multithreaded compilation if OMP is set above
## e.g.(-openmp for intel, -fopenmp for gnu)
#FFLAGS_OMP

## Other compiler flags
#FFLAGS_OTHER

## Override name of linker (default is to use fortran compiler)
#LD

##Name of c compiler (for svnrev version tracker)
#CC = cc

## Name for GKW if the version number cannot be found with svnrev
#DEFAULT_GKW_VERSION = SOME_VERSION
