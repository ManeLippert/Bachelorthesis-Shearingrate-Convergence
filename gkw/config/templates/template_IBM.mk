## --  Template for IBM Fortan compiler.  See also ../babel and ../vip

#REQUIREDMODULES = IBM Fortran compiler, FFTW, MPI

MPI = usempi2
FFTLIB =FFT_FFTW3

CC = cc
FC = mpxlf95_r  # or mpxlf90_r

FFLAGS_OPT = -O4
#FFLAGS_OPT = -q64 -O2 -qnoipa  -qmaxmem=-1 -qflag=I:I

FFLAGS_DOUBLE =-qautodbl=dbl4 # can try also =dbl

# DEBUG=on
# If DEBUG=on, might want only a subset of these
FFLAGS_DEBUG = -g -qsigtrap -qoptdebug -qfullpath -C -qinitauto=7FBFFFFF \
      -qflttrap=overflow:underflow:zerodivide:invalid:enable \
      -qfloat=nans -qsigtrap -qkeepparm

#SMP = OPENMP
FFLAGS_OMP = -qsmp=omp

# IBM fortran compiler has unusual preprocessing syntax
PREPROC_FLAG =-WF
PREPROC_SEP = ,
PREPROC_PREFIX =-D

# Compiling umfpack library can require different C preprocessor flags
# See ../turing for a full example
#C_PREPROC_FLAG= 
#C_PREPROC_SEP=${space}
#C_PREPROC_PREFIX = -D

FFLAGS_INC = -I. -I${FFTW_HOME}/include
LDFLAGS = -L${FFTW_HOME}/lib -lfftw3 -lfftw3f

# Make code timings, no reason not to use
PERF = perf 

# Need the lines below if you want to use nonspectral or implicit scheme
# GKW_LIBS = /path/to/compiled/gkw/libs   # compile them with the same compiler 
# IMPLICIT=umfpack
# LDFLAGS+=-L${GKW_LIBS}/UMFPACK/Lib/ -lumfpack_gkw64_NB -L${GKW_LIBS}/AMD/Lib/ -lamd

