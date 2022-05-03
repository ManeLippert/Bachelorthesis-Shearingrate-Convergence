## Name of fortran compiler (default is gfortran)
FC = /home/bt161087/mpich/bin/mpif90 
## Compiler debugging flags used if DEBUG=on (e.g. -g -C -traceback)
FFLAGS_DEBUG = -g -fbounds-check -fbacktrace -Wunused #-ffpe-trap=zero,overflow #-Wall -W -Wunderflow

## Compiler optimisation flasgs used if OPTFLAGS=on (e.g -03)
FFLAGS_OPT = -O3 -ftracer -fomit-frame-pointer -pipe -fweb

## Compiler flag used when compiling in double precision (default is for gfortran)
FFLAGS_DOUBLE = -fdefault-real-8 -fdefault-double-8

## Compiler include flags ( e.g. -I/path/to/include)
FFLAGS_INC = -I. -I/home/bt161087/mpich/include -I/usr/include 
## Flags required for multithreaded compilation if SMP is set above
## e.g.(-openmp for intel, -fopenmp for gnu)
FFLAGS_OMP= -fopenmp
## Other compiler flags
FFLAGS_OTHER = -std=f95
