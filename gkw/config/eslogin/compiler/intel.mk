# -- Archer HPC UK resource: Intel compiler

REQUIREDMODULES += PrgEnv-intel 
FC = ftn
#FFLAGS_WARN  = -std95
FFLAGS_OPT = -O3 #-ipo
FFLAGS_DOUBLE = -r8
FFLAGS_OMP = -openmp
FFLAGS_DEBUG = -g -C -traceback
