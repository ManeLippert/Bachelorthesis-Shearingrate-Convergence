## -- Archer HPC UK resource: Cray compiler

REQUIREDMODULES += PrgEnv-cray 
FC = ftn
FFLAGS_OPT = -O3 
FFLAGS_DOUBLE = -s real64 # -default64 fails in mpi_cart_coords
FFLAGS_OMP = -h omp #enabled in ftn wrapper by default.  To disable, need -h noomp
FFLAGS_DEBUG = -g
