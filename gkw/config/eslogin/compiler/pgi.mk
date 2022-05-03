## -- Archer HPC UK resource: Portland compiler --

REQUIREDMODULES = PrgEnv-pgi
CC = gcc
FC = ftn
FFLAGS_OPT = -fastsse 
FFLAGS_DOUBLE = -r8
FFLAGS_DEBUG = -traceback -gopt
FFLAGS_OMP = -mp=nonuma -Minfo=mp
