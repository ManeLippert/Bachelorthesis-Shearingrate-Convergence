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
OPTFLAGS = off

# Use the gnu compiler settings in subfolder ass the default
COMPILER=gnu

## Include the HDF5 libraries to write in the HDF5 Format.
IO_LIB = HAVE_HDF5

## Include the slepc and petsc libraries for the eigenvalue solver.
## If the libraries are not in the search path, you must also set
## LDFLAGS and FFLAGS_INC appropriately
SLEPC = HAVE_SLEPC

## Include the librsb library with (potentially faster)
## sparse linear algebra routines
LIBRSB = HAVE_LIBRSB

## Include the MKL library with (potentially faster)
## sparse linear algebra routines
#LIBMKL = HAVE_MKL

## Those flags were copied from $(h5fc -show)
FFLAGS_H5 = -L/home/btpp/btp00000/shared/hdf5-current/lib64 \
            -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 \
            -L/home/btpp/btp00000/shared/hdf5-current/lib -lsz -lz -ldl -lm \
            -Wl,-rpath -Wl,/home/btpp/btp00000/shared/hdf5-current/lib64
FFLAGS_H5_INC = -I/home/btpp/btp00000/shared/hdf5-current/include

# Note that the SLEPC/PETC library depends also on MPI
FFLAGS_SLEPC = -lslepc -lpetsc -L/home/btpp/btp00000/shared/petsc-current/lib -L/home/btpp/btp00000/shared/slepc-current/lib
FFLAGS_SLEPC_INC = -I/home/btpp/btp00000/shared/petsc-current/include -I/home/btpp/btp00000/shared/slepc-current/include

FFLAGS_FFTW3 =  -lfftw3 -lfftw3f -L/usr/lib64 
FFLAGS_FFTW3_INC = -I/usr/include

FFLAGS_MPI = -lmpichf90 -L/home/btpp/btp00000/shared/mpich-current/lib64 -lmpifort -Wl,-rpath -Wl,/home/btpp/btp00000/shared/mpich-current/lib64 -Wl,--enable-new-dtags -lmpi -lrt -lpthread # -lmpich -lrt -lopa -lmpl
FFLAGS_MPI_INC = -I/home/btpp/btp00000/shared/mpich-current/include

FFLAGS_UMFPACK = -L../../libs/UMFPACK/Lib/ -lumfpack_gkw64_NB -L../../libs/AMD/Lib/ -lamd
FFLAGS_UMFPACK_INC = #none

# Directory where the MKL library is located
# FFLAGS_MKL_INC= -I$(MKL_HOME)/include/intel64/lp64  -I${MKL_HOME}/include
# FFLAGS_MKL = -L "${MKL_HOME}/lib/intel64" -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -L "${MKL_HOME}/../compiler/lib/intel64" -liomp5 -lm
# MKL_SPBLAS_MODULE_SRC = ${MKL_HOME}/include/mkl_spblas.f90


## Those flags were copied from $(librsb-config --ldflags --extra_libs)
FFLAGS_LIBRSB =  -L/home/btpp/btp00000/shared/librsb-current/lib64 -Wl,-rpath -Wl,/home/btpp/btp00000/shared/librsb-current/lib64 -lrsb -lm -lgfortran
#FFLAGS_RSB = $(shell librsb-config --I_opts --static --ldflags --extra_libs)
# -L/home/btpp/btp00000/shared/librsb-current/lib64 
FFLAGS_LIBRSB_INC = -I/home/btpp/btp00000/shared/librsb-current/include

FFLAGS_OMP = -fopenmp

## The complete set of include flags
FFLAGS_INC = -I. $(FFLAGS_H5_INC) $(FFLAGS_MPI_INC) $(FFLAGS_SLEPC_INC) $(FFLAGS_UMFPACK_INC) $(FFLAGS_FFTW3_INC) $(FFLAGS_LIBRSB_INC) $(FFLAGS_MKL_INC)

## The complete set of linking flags ( e.g. -L/path/to/include -lfftw3)
LDFLAGS =  $(FFLAGS_MPI) $(FFLAGS_SLEPC) $(FFLAGS_H5) $(FFLAGS_FFTW3) $(FFLAGS_UMFPACK) $(FFLAGS_LIBRSB) $(FFLAGS_MKL)


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
