## ---------- RZG-VIP-Power6--------------

MPI = mpi2
FFTLIB =FFT_FFTW3

CC = cc
FC = mpxlf95_r
FFLAGS_OPT = -q64 -O2 -qnoipa  -qmaxmem=-1 -qflag=I:I
FFLAGS_DOUBLE =-qautodbl=dbl4
PREPROC_FLAG =-WF
PREPROC_SEP = ,
PREPROC_PREFIX =-D
FFLAGS_INC = -I. -I/u/system/Power6/libs/fftw/fftw-3.1.2/include/
LDFLAGS = -L/u/system/Power6/libs/fftw/fftw-3.1.2/lib -lfftw3 -lfftw3f
