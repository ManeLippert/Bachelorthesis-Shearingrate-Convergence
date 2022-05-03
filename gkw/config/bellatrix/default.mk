## -- Warwick desktop AGP
MPI = mpi2
FFTLIB = FFT_FFTW3
IMPLICIT = umfpack
#SMP = OPENMP

FC=mpif90
FFLAGS_OPT= -O2 -no-prec-div 
FFLAGS_DOUBLE = -r8
FFLAGS_OTHER = -vec-report0
LDFLAGS = -L/home/space/phsfar/effe/gkw/UMFPACK/Lib -lsvml /home/space/phsfar/effe/gkw/UMFPACK/Demo/umf4_f77wrapper64.o -lumfpack -lamd -lm -lfftw3 -lfftw3f
