# GKW basic makefile - for use when the more comprehensive GNU makefile
# cannot be used.
#
#============================================================================
### EXECUTABLE name
### (change the name here if you wish)
EXEC=../gkw.x

#============================================================================
### PREPROCESSOR OPTIONS
### (it is unlikely you need to modify anything beyond the first 2 options
###  below)

## --- FFT library ---
## The prefered FFT library is FFTW3; interfaces for some others may
## exist in the fft module of the code. To switch off FFT's altogther
## via a dummy interface, set `FFTLIB=nofft'.
FFTLIB=FFT_FFTW3
#FFTLIB=nofft

## --- Parallel switch ---
## For compiling without MPI, select `PARALLEL=nompi'.
PARALLEL=mpi2
#PARALLEL=nompi

## --- Set the *real* precision ---
## This should probably remain unchanged. If you really wish to run in
## single precision, change this to `real_precsion_default', then adjust
## the compiler precision flags in FFLAGS below. In addition, you may have
## link with FFTW3 differently if you are using that option for FFTLIB above
## (usually -lfftw3f, rather than -lfftw3). It appears that GKW will run
## successfully in single precision, therefore saving some time and memory.
## However, we presently give no assurances that this will work.
REAL_PRECISION=real_precision_double
#REAL_PRECISION=real_precision_default

## --- Version switch ---
## Do not need to change this.
GKW_VERSION = NOVERSION

## NO NEED TO MODIFY THE LINE BELOW UNLESS MORE OPTIONS ARE ADDED
## (slepc, petsc, umfpack, mpif90_interface) not listed here
FPREPROCS = -D$(FFTLIB)  \
			-D$(PARALLEL) \
			-D$(REAL_PRECISION) \
			-D$(GKW_VERSION) \

#============================================================================
### COMPILER OPTIONS

## --- FORTRAN compiler ---
## The default here is a wrapper provided by most MPI implementations.
#FC=mpif90
FC=gfortran

## --- FORTRAN compiler flags ---
## We use double precision as the default, so please replace `-r8' with the
## appropriate switch for your compiler (gfortran uses -fdefault-real-8).
## Most compilers support `-O2' to turn on reasonable optimisation.
FFLAGS= -fdefault-real-8 -O2 -DUSE_MPIF -fallow-invalid-boz -fallow-argument-mismatch -fno-range-check
#FFLAGS= -fdefault-real-8 -O2 #gfortran

## --- Include PATH ---
## You may need to add something here if the MPI, FFTW or other library header
## files are installed in an unusual location.
#INC=
INC= -I/usr/include -I/usr/lib/x86_64-linux-gnu/fortran/gfortran-mod-15/openmpi
#INC= -I/usr/local/some_directory/include/

## --- LINKER ---
## (usually, the FORTRAN compiler should be fine for linking)
LD=$(FC) $(FFLAGS)

## --- linking flags ---
## It may be necessary to add flags here if you did not use the `mpif90',
## `h5pfc' or similar wrapper for FC. If you do not want to use FFTW3 you
## can (and should) remove the default flag here. As with the include PATH,
## if your library is in an unusual place, you may need to set that path here.
LDFLAGS = -lfftw3 -lfftw3f
#LDFLAGS= -L/usr/local -lfftw3
#LDFLAGS= -L/usr/local/some_directory/lib -lfftw3f

#############################################################################
#############################################################################
#---------- YOU SHOULD NOT NEED TO EDIT ANYTHING BEYOND HERE ---------------#
#           (unless more sources are added to the code)                     #
#############################################################################
#############################################################################

#============================================================================
### FORTRAN object list

include objfiles.mk

#============================================================================
### rules

.SUFFIXES:
.SUFFIXES: .f90 .F90 .o

.f90.o:
	    $(FC) $(FFLAGS) $(INC) -o $*.o -c $*.f90
.F90.o:
	    $(FC) $(FFLAGS) $(INC) $(FPREPROCS) -o $*.o -c $*.F90

#============================================================================
### various targets

default: $(EXEC)

$(EXEC) : $(OBJLIST)
	$(LD) $(LDFLAGS) $(OBJLIST) -o $(EXEC)

clean:
	  rm -f *.o external/*.o *.mod work.pc work.cc1 gkw_info.h $(EXEC)

gkw_info.h:

#============================================================================
### required FORTRAN deps

include deps.mk
