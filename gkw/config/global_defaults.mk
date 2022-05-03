## This file is always included by the gkw makefile. The variables set here
## are the defaults which are overridden by the values in the
## host/user/compiler specific file(s). See template.mk for details of the
## variables.

OPTFLAGS = on
DEBUG    = off
WARN     = on

REAL_PRECISION = real_precision_double
#REAL_PRECISION = real_precision_default

CC     = cc
CFLAGS = -O2
FC = gfortran
FFLAGS_DOUBLE = -fdefault-real-8 -fdefault-double-8
FFLAGS_OPT = -O2
FFLAGS_DEBUG = -g
FFLAGS_INC = -I.
#FFLAGS_WARN = -W -Wall

UMFPACK_INTERNAL=$(LIBDIR)/UMFPACK/Lib/libumfpack.a $(LIBDIR)/AMD/Lib/libamd.a
LD_UMFPACK=-L$(LIBDIR)/UMFPACK/Lib -lumfpack -L$(LIBDIR)/AMD/Lib -lamd

INPUT_CHECK = input_check
DEFAULT_GKW_VERSION = UNKNOWN

RANLIB=touch

#APS: some config files need GKW_LIBS
GKW_LIBS=$(GKW_HOME)/libs
