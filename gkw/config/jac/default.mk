## -- JET: jac 32 bit machines- gfortran --
## with hd5f and umfpack

MPI = usempi2
FFTLIB = FFT_FFTW3
SMP = OPENMP
#DEBUG = on

# you have to load gfortran first !
REQUIREDMODULES=standard gfortran openmpi hdf5

FC = mpif90
FFLAGS_INC = -I. -I/usr/include -I${FFTW_HOME}/include
FFLAGS_OTHER = -std=f95 #-W -Wall -Wunderflow
FFLAGS_OMP = -fopenmp
FFLAGS_OPT = -O3 -ftracer -fomit-frame-pointer -pipe -fweb
FFLAGS_DEBUG = -g -fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow
FFLAGS_DOUBLE = -fdefault-real-8 -fdefault-double-8
LDFLAGS = -lfftw3 -lfftw3f

# Need the lines below if you want to use nonspectral, global or implicit scheme
IMPLICIT=umfpack32
LDFLAGS+=-L${GKW_LIBS}/UMFPACK/Lib/ -lumfpack_gkw32_NB -L${GKW_LIBS}/AMD/Lib/ -lamd

# Need the lines below if you want to use hdf5 output format options
IO_LIB = HAVE_HDF5
HDF5_HOME=/usr/local/depot/hdf5-1.8.11-gfortran
LDFLAGS+= -L/${HDF5_HOME}/lib -lhdf5 -lhdf5_fortran -lhdf5_hl -lhdf5hl_fortran -lz -lm
FFLAGS_INC+= -I/${HDF5_HOME}/include
