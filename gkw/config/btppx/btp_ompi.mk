# Bayreuth desktop generic - gfortran and homebuilt openmpi-1.4.3
# To use with the correct mpirun, do
# export PATH="/home/btpp/btpp04/openmpi-1.4.3/bin:$PATH"
# To get the correct libraries at runtime, you also need
# export LD_LIBRARY_PATH=''
 
## Without MPI, only single processor runs can be performed 
## Uncomment to use MPI:
MPI = mpi2

## Without FFTW (default), only linear runs can be performed (off by default) 
## Uncomment to use FFTW:
FFTLIB = FFT_FFTW3

## Turn compiler debugging on (off by default)
#DEBUG = on

## Turn compiler optismations off (on by default)
#OPTFLAGS = off

## Name of fortran compiler (default is gfortran)
# RB: Could not get this working, due to 'error while loading shared
#     libraries: libopen-pal.so.0: cannot open shared object file: No such
#     file or directory'.
#     (Changing LD_LIBRARY_PATH also does not work, other files can not be
#     found either).
#     The whole directory seems to have been moved. (From /home/btpp/btpp04/.)
#FC = /home/btpp/bt482456/btpp04/openmpi-1.4.3/bin/mpif90
FC = mpif90

## Compiler debugging flags used if DEBUG=on (e.g. -g -C -traceback)
FFLAGS_DEBUG = -g -fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow -Wall -W -Wunderflow

## Compiler optimisation flasgs used if OPTFLAGS=on (e.g -03)
FFLAGS_OPT = -O3 -ftracer -fomit-frame-pointer -pipe -fweb

## Compiler flag used when compiling in double precision (default is for gfortran)
FFLAGS_DOUBLE = -fdefault-real-8 -fdefault-double-8

## Compiler include flags ( e.g. -I/path/to/include)
FFLAGS_INC = -I. -I/home/btpp/bt482456/btpp04/openmpi-1.4.3/include -I/usr/include 

## Linking flags ( e.g. -L/path/to/include -lfftw3)
LDFLAGS = -L/usr/lib64 -lfftw3 -lfftw3f -L/home/btpp/bt482456/btpp04/openmpi-1.4.3/lib/ -Wl,-rpath,/home/btpp/bt482456/btpp04/openmpi-1.4.3/lib/

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

## Flags required for multithreaded compilation if SMP is set above
## e.g.(-openmp for intel, -fopenmp for gnu)
FFLAGS_OMP= -fopenmp

## Other compiler flags
#FFLAGS_OTHER

## Override name of linker (default is to use fortran compiler)
#LD

##Name of c compiler (for svnrev version tracker)
#CC = cc

## Name for GKW if the version number cannot be found with svnrev
#DEFAULT_GKW_VERSION = SOME_VERSION
