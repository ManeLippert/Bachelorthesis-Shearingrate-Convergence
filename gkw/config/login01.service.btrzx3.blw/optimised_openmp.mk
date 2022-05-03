## --- Btcluster with out floating point exception checks -> compiler optimisations ---
# Inherit defaults

# by default, some runtime checks are enabled, to force more people to
# notice errors without convincing them to use certain make parameters
include ${CONFIGDIR}/default.mk

SMP = OPENMP
FFLAGS_OMP = -qopenmp #-parallel  -xAVX -Wall -vec-report3 -opt-report3 -restrict -guide

# for faster production runs on btcluster these runtime checks can be disabled.
OPTFLAGS = on
DEBUG = off
FFLAGS_OTHER =
# In order to do this, compile with
# gkwmake -j CONFIG=config/login01.service/optimised_openmp.mk DEBUG=off OPTFLAGS=on
