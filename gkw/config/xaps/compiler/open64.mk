# --- Open64 fortran compiler on xaps ---

FC              = mpif90-open64
FFLAGS_INC      = -I. -I/usr/include
FFLAGS_OTHER    = 
FFLAGS_OPT      = -O3
FFLAGS_OMP      = -openmp
FFLAGS_DEBUG    = -g
FFLAGS_DOUBLE   = -r8

