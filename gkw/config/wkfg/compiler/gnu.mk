# -- gfortran with openMPI wrapper --
FC = mpif90
FFLAGS_DEBUG = -g -fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow
FFLAGS_OPT = -O3 -ftracer -fweb -pipe
FFLAGS_DOUBLE = -fdefault-real-8 -fdefault-double-8
FFLAGS_OMP = -fopenmp
FFLAGS_INC = -I. -I/usr/include
FFLAGS_OTHER = -march=native -W -Wall -Wunderflow
LD = mpif90
