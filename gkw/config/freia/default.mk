## -- JET: jac 64 bit machines - gfortran --
## with hd5f and umfpack

MPI = usempi2
FFTLIB = FFT_FFTW3
SMP = OPENMP
#DEBUG = on

# you have to load gfortran then hdf5 first !
REQUIREDMODULES=gfortran hdf5 openmpi

FC = mpif90
FFLAGS_INC = -I. -I/usr/include #-I${FFTW_HOME}/include
FFLAGS_OTHER = #-std=f95 #-W -Wall -Wunderflow
FFLAGS_OMP = -fopenmp
FFLAGS_OPT = -O3 -ftracer -fomit-frame-pointer -pipe -fweb
FFLAGS_DEBUG = -g -fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow
FFLAGS_DOUBLE = -fdefault-real-8 -fdefault-double-8
LDFLAGS = -lfftw3 -lfftw3f

# Need the lines below if you want to use nonspectral, global or implicit scheme
IMPLICIT=umfpack
#GKW_LIBS=/home/fcasson/gkw/libs64
LD_UMFPACK=-L${GKW_LIBS}/UMFPACK/Lib/ -lumfpack -L${GKW_LIBS}/AMD/Lib/ -lamd

# Need the lines below if you want to use hdf5 output format options
#IO_LIB = HAVE_HDF5
#LDFLAGS+= -lhdf5 -lhdf5_fortran -lhdf5_hl -lhdf5hl_fortran -lz -lm

