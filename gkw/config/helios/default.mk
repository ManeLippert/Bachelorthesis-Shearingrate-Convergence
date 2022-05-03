 ## --- IFERC helios CSC ---
REQUIREDMODULES = intel bullxmpi fftw/3.3
MPI = usempi2
FFTLIB = FFT_FFTW3
SMP = OPENMP

FC = mpif90
#To use intelmpi instead, first load module for intelmpi instead of bullxmpi
#and then export USE_IMPI=1 for gkwnlin
#FC = mpiifort

FFLAGS_DEBUG = -fpe0 -O0 -g -traceback -check assume,pointers,stack,uninit,bounds
FFLAGS_OPT = -axCORE-AVX2 -O3 -no-prec-div -ip # slower compile and no speedup
FFLAGS_OPT = -O2 -ip #-no-prec-div #-xAVX #problems in umfpack, and no faster
FFLAGS_OMP = -openmp
FFLAGS_DOUBLE = -r8
FFLAGS_INC = -I/${FFTW_DIR}/include
LDFLAGS = -L/${FFTW_DIR}/lib -lfftw3 #-lfftw3f

#Comment if umfpack is not required
#GKW_LIBS = /csc/home1/fjc/gkw/libs # shaerd permissions on home no longer allowed
#GKW_LIBS =  /csc/workdir1/fjc/gkw_libs # copied to work, should be readable
#LDFLAGS+= -L${GKW_LIBS}/UMFPACK/Lib/ -lumfpack_gkw64_NB -L${GKW_LIBS}/AMD/Lib/ -lamd 

# or simple now to build your own version of umfpack
# linking is also simpler
IMPLICIT=umfpack
LDFLAGS+= -lumfpack

