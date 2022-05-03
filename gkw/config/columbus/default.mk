## -- CCFE CULHAM cluster (very old software, not updated since 2009)  --

# There is a version of mpich2 installed, which only seems
# to work with gfortran.  The gfortran libriares do not 
# appear in the path (there seems to be no module to load 
# GNU compiler path settings), so some tricks with softlinks
# are used instead (but may need updating in future).


MPI = mpi2
FFTLIB = FFT_FFTW3
REQUIREDMODULES = mpich2 
CC = gcc
FC = mpif90
FFLAGS_OPT = -O2
FFLAGS_DEBUG = -g  
#FFLAGS_OMP=-fopenmp
FFLAGS_DOUBLE = -fdefault-real-8
FFLAGS_INC = -I/usr/local/fusion/64/include -I./
LDFLAGS = -L/usr/local/fusion/64/lib -lfftw3 -L/lib64 -L../../config/columbus -Wl,-rpath,${GKW_HOME}/config/columbus
#SMP=OPENMP
#IMPLICIT=umfpack
