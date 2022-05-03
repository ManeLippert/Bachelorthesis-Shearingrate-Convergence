## Name of fortran compiler (default is gfortran)
FC = gfortran
## Compiler debugging flags used if DEBUG=on (e.g. -g -C -traceback)
FFLAGS_DEBUG = -g -Og -fbacktrace  #-Wall

## Compiler optimisation flasgs used if OPTFLAGS=on (e.g -03)
FFLAGS_OPT = -O2 -ftracer -fomit-frame-pointer -pipe -fweb -march=native

## Compiler flag used when compiling in double precision (default is for gfortran)
FFLAGS_DOUBLE = -fdefault-real-8 -fdefault-double-8

## Flags required for multithreaded compilation if SMP is set above
## e.g.(-openmp for intel, -fopenmp for gnu)
FFLAGS_OMP= -fopenmp 
## Other compiler flags
FFLAGS_OTHER = -std=f2008 -Wunused -Wcompare-reals -ffpe-trap=zero,overflow,invalid -fcheck=bounds,do,mem,pointer,recursion  -fmax-identifier-length=31 -Wno-tabs -DFORCE_IO_ATOMICITY
# -Wcompare-reals can be activated with gfortran 3.8 or newer, not available in 3.7.2.

#FFLAGS_WARN = -W -Wall
