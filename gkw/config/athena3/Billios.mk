## Without MPI, only single processor runs can be performed 
## Uncomment to use MPI:
MPI = mpi2

## Without FFTW (default), only linear runs can be performed (off by default) 
## Uncomment to use FFTW:
FFTLIB = FFT_FFTW3

## Turn compiler debugging on (off by default)
DEBUG = on

## Turn compiler optismations off (on by default)
OPTFLAGS = off

# Use the gnu compiler settings in subfolder ass the default
#COMPILER=gfortran

## Include the HDF5 libraries to write in the HDF5 Format.
#IO_LIB = HAVE_HDF5

## Include the slepc and petsc libraries for the eigenvalue solver.
## If the libraries are not in the search path, you must also set
## LDFLAGS and FFLAGS_INC appropriately
#SLEPC = HAVE_SLEPC

# Note that the SLEPC/PETC library depends also on MPI
#FFLAGS_SLEPC = -lslepc -lpetsc -L/home/btpp/btp00000/shared/petsc-current/lib -L/home/btpp/btp00000/shared/slepc-current/lib
#FFLAGS_SLEPC_INC = -I/home/btpp/btp00000/shared/petsc-current/include -I/home/btpp/btp00000/shared/slepc-current/include

FFLAGS_FFTW3 = -L/usr/local/lib -lfftw3 #-lfftw3f 
FFLAGS_FFTW3_INC = -I/usr/local/include

# tradionally used settings:
FFLAGS_MPI = -L/usr/local/lib -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi 
FFLAGS_MPI_INC = -I/usr/local/include -Wl,-flat_namespace -Wl,-commons,use_dylibs

FFLAGS_UMFPACK = -L../../libs/UMFPACK/Lib/ -lumfpack_gkw64_NB -L../../libs/AMD/Lib/ -lamd
FFLAGS_UMFPACK_INC = #none

## The complete set of include flags
#FFLAGS_INC =  $(FFLAGS_H5_INC) $(FFLAGS_MPI_INC) $(FFLAGS_SLEPC_INC) $(FFLAGS_UMFPACK_INC) $(FFLAGS_FFTW3_INC) -I.
FFLAGS_INC = -I. -I/usr/local/include

## The complete set of linking flags ( e.g. -L/path/to/include -lfftw3)
LDFLAGS =  $(FFLAGS_MPI) $(FFLAGS_SLEPC) $(FFLAGS_H5) $(FFLAGS_FFTW3)  $(FFLAGS_UMFPACK) -lgomp

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
#PERF = perf  # performance timings
#PREPROC_FLAG
#PREPROC_SEP
#PREPROC_PREFIX

## Override name of linker (default is to use fortran compiler)
#LD

##Name of c compiler (for svnrev version tracker)
#CC = cc

## Name for GKW if the version number cannot be found with svnrev
#DEFAULT_GKW_VERSION = SOME_VERSION
