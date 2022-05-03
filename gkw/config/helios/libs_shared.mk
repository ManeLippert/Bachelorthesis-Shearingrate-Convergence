 ## --- IFERC helios CSC ---
REQUIREDMODULES = intel bullxmpi fftw/3.3
MPI = usempi2
FFTLIB = FFT_FFTW3
SMP = OPENMP

FC = mpif90
#To use intelmpi instead, first load module for intelmpi instead of bullxmpi
#and then export USE_IMPI=1 for gkwnlin
#FC = mpiifort

FFLAGS_DEBUG = -fpe0 -C -g -traceback
FFLAGS_OPT = -axCORE-AVX2 -O3 -no-prec-div -ip # slower compile and no speedup
FFLAGS_OPT = -O2 -ip -no-prec-div #-xAVX #problems in umfpack, and no faster
FFLAGS_OMP = -openmp
FFLAGS_DOUBLE = -r8
FFLAGS_INC = -I/${FFTW_DIR}/include
LDFLAGS = -L/${FFTW_DIR}/lib -lfftw3 #-lfftw3f

#use precompiled UMFPACK from $PROJECT dir
GKW_LIBS=/csc/project/RESIDUAL/libs# should be readable to members
UMFPACK_LIBRARY = pre-built
IMPLICIT=umfpack
LD_UMFPACK = -lumfpack -lamd -L${GKW_LIBS}/AMD/Lib -L${GKW_LIBS}/UMFPACK/Lib
