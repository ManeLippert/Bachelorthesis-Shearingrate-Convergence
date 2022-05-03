# Bayreuth desktop generic - gfortran and homebuilt mpich2
 
## Without MPI, only single processor runs can be performed 
## Uncomment to use MPI:
MPI = usempi2

## Without FFTW (default), only linear runs can be performed (off by default) 
## Uncomment to use FFTW:
FFTLIB = FFT_FFTW3

## Turn compiler debugging on (off by default)
DEBUG = on

## Turn compiler optismations off (on by default)
#OPTFLAGS = off

# Use the gnu compiler settings in subfolder ass the default
COMPILER=gnu

## Linking flags ( e.g. -L/path/to/include -lfftw3)
LDFLAGS = -L/usr/lib64 -lfftw3 -lfftw3f  -L/home/bt161087/mpich/lib64 -lmpich  -lrt -lpthread -L../../libs/UMFPACK/Lib/ -lumfpack_gkw64_NB -L../../libs/AMD/Lib/ -lamd

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
IMPLICIT = umfpack
#HDF5 = HDF5
PERF = perf  # performance timings
#PREPROC_FLAG
#PREPROC_SEP
#PREPROC_PREFIX

## Override name of linker (default is to use fortran compiler)
#LD

##Name of c compiler (for svnrev version tracker)
#CC = cc

## Name for GKW if the version number cannot be found with svnrev
#DEFAULT_GKW_VERSION = SOME_VERSION
