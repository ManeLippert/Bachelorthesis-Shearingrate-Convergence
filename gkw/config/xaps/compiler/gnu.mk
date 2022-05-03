# --- GNU fortran compiler on xaps ---

REQUIREDMODULES_IS_BROKEN = mpich2-1.0.7-gnu-4.3.1-x86_64

FC              = mpif90
FFLAGS_INC      = -I. -I/usr/include
FFLAGS_OTHER    = -W -Wall -Wunderflow
FFLAGS_OMP      = -fopenmp
FFLAGS_OPT      = -O3 -ftracer -fomit-frame-pointer -pipe -fweb
FFLAGS_DEBUG    = -g -fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow
FFLAGS_DOUBLE   = -fdefault-real-8 -fdefault-double-8
