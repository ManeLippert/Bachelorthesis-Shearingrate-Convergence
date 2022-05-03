# --- ITER cluster ----
# 
# Need to unload UCX and numactl modules to avoid UCX error for parallel runs
#
# Test cases ok, except:
#     global_em in serial (SEG fault, not sure why)
# 

REQUIREDMODULES= FFTW/3.3.8-intel-2020a # and then unload UCX and numactl


MPI = usempi2
FFTLIB = FFT_FFTW3
#SMP = OPENMP

OPTFLAGS = on
DEBUG = off
PERF=perf

CC = icc
FC = mpiifort
#CC = gcc
#FC = mpifort
LD = $(FC)

IMPLICIT=umfpack
#IO_LIB = HAVE_HDF5
GKW_LIBS=$GKW_HOME/libs

# intel 
FFLAGS_DOUBLE = -r8
FFLAGS_OPT = -O3 -ip -noalign
FFLAGS_OMP = -qopenmp
FFLAGS_DEBUG = -fpe0 -O0 -g -traceback -check assume,pointers,stack,uninit,bounds
# gfortran
#FFLAGS_DOUBLE = -fdefault-real-8 -fdefault-double-8
#FFLAGS_OPT = -O3 -ftracer -fomit-frame-pointer -pipe -fweb
#FFLAGS_OMP = -fopenmp
#FFLAGS_DEBUG = -g -fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow

#FFLAGS_INC =`pkg-config fftw3 hdf5 PETSc SLEPc zlib --cflags`
#FFLAGS_INC = -I. `pkg-config fftw3 --cflags`
FFLAGS_INC = -I. -I/work/imas/opt/EasyBuild/software/FFTW/3.3.8-intel-2020b/include
#FFLAGS_INC = -I.
#LDFLAGS = `pkg-config fftw3 hdf5 PETSc SLEPc zlib --libs`
LDFLAGS = -L/work/imas/opt/EasyBuild/software/FFTW/3.3.8-intel-2020b/lib -lfftw3 -lfftw3f
#LDFLAGS = -lfftw3 


