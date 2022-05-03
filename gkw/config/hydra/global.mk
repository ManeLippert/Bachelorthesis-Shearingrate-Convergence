## --- IPP hydra with umfpack implicit library  ---
# Inherit defaults
include ${CONFIGDIR}/default.mk

#GKW_LIBS=/u/wiho/gkw_git2/libs
IMPLICIT=umfpack
LDFLAGS+=-L/$(GKW_LIBS)/UMFPACK/Lib/ -lumfpack_gkw64_NB -L/${GKW_LIBS}/AMD/Lib/ -lamd
