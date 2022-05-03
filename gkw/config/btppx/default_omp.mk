# Bayreuth desktop generic - gfortran and homebuilt mpich2 + openmp support.

# start with default setting.
include ${CONFIGDIR}/default.mk

# Now the changes to the defaults:
# Run without debug information and with optimization
DEBUG = off
OPTFLAGS = on

# Enable also openmp.
SMP = OPENMP
