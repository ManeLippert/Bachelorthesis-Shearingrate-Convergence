##############################################################################
### Makefile for the BlueGene/Q Turing @ IDRIS
##############################################################################
# use: module load fftw/3.3.2

OPTFLAGS = on
DEBUG = off

# IBM compilers
FC = mpixlf95_r
F77=bgxlf
CC=bgcc
RANLIB=ranlib

FFLAGS_OPT = -O4
FFLAGS_DEBUG = -g -qsigtrap -qoptdebug -qfullpath -C -qinitauto=7FBFFFFF \
      -qflttrap=overflow:underflow:zerodivide:invalid:enable \
      -qfloat=nans -qsigtrap -qkeepparm 

FFLAGS_DOUBLE = -qautodbl=dbl4

FFLAGS_INC = -I.

FFLAGS_OMP = -qsmp=omp

#FFLAGS_OTHER = -qarch=450d

# Umfpack compiles correctly with the new gkw procedure, the compile_umf script is not needed
# IMPLICIT=umfpack
# UMFPACK_LIBRARY=compile
#LD_UMFPACK=-lumfpack -L${GKW_HOME}/libs/UMFPACK/Lib
#It appears to be necessary to explicitly link -lamd (not sure why)
#LD_UMFPACK = -L${GKW_HOME}/libs/UMFPACK/Lib -lumfpack -L${GKW_HOME}/libs/AMD/Lib -lamd

LD = $(FC)

# The IBM Fortran compiler requires strange preprocessor options
PREPROC_FLAG   = -WF
PREPROC_SEP    = ,
PREPROC_PREFIX = -D

# But, for umfpack compilation, the C preprocessing should be set as usual (at least with gcc)
# note, "space" can only be specified with the variable, since whitespace is ignored
C_PREPROC_FLAG= 
C_PREPROC_SEP=${space}
C_PREPROC_PREFIX = -D

REAL_PRECISION = real_precision_double

MPI = mpi2

FFTLIB = FFT_FFTW3

SMP = OPENMP

PERF = perf2 

MPIRUNCMD = mpirun
