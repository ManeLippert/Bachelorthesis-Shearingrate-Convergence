## ---------- JET jac* - pgi --------------
##
MPI = mpi2
FFTLIB = FFT_FFTW3
REQUIREDMODULES=standard pgi openmpi
CC = cc
FC = mpif90
FFLAGS_OPT = -fastsse
FFLAGS_DOUBLE =-r8
FFLAGS_OMP=-mp
FFLAGS_DEBUG=-g -C
LDFLAGS = -lfftw3 -lfftw3f
#SMP=OPENMP
#PERF=perf # non neglible overhead
