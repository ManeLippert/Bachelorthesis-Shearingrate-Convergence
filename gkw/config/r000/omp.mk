## --- CINECA Marconi with OpenMP shared memory parallelism ---
# Inherit defaults
include ${CONFIGDIR}/default.mk

# These are for OpenMP
SMP = OPENMP
FFLAGS_OMP = -fopenmp
