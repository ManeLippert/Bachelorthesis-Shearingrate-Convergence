##############################################################################
###
###    GKW src GNU makefile
###
###
##############################################################################

## Clear out any defaults
.SUFFIXES:
.DEFAULT:
LD=

## where to look for various *source* files
vpath          # empty to start with; add more explicitly as needed
vpath %.F90    $(SRCDIR)
vpath %.f90    $(SRCDIR)
vpath %.c      $(SRCDIR)

## Include the defaults for switches, compilers, flags and preprocessing.
include $(DEFAULTS)

## The specific configuration file which has been selected via hostname,
## username etc. or otherwise set explicitly when running make.
include $(CONFIGFILE) $(COMPILERFILE)

## (some variables used below)
empty:=
space:= $(empty) $(empty)
colon:= $(empty):$(empty)

##############################################################################
### Preprocessor option manipulation (see $(TEMPLATE) as set in the top level#
### makefile for used variables). Variables are filtered here in the makefile#
### before being being passed to the compiler. In some cases, blank or dummy #
### values are used when an unknown option is provided; these can be         #
### converted to errors where necessary (e.g. in the case of REAL_PRECISION, #
### which is always required).                                               #
##############################################################################


#### REAL_PRECISION -> _REAL_PRECISION
# allowed input values:   "real_precision_default", "real_precision_double"
# preprocessor variables: "real_precision_default", "real_precision_double"
ifneq ($(_REAL_PRECISION),)
  $(error _REAL_PRECISION was passed to make or set via an environment variable)
else
  ifneq ($(REAL_PRECISION),)
    ifeq ($(REAL_PRECISION),real_precision_default)
      _REAL_PRECISION = $(REAL_PRECISION)
    else
      ifeq ($(REAL_PRECISION),real_precision_double)
        _REAL_PRECISION = $(REAL_PRECISION)
      else
        $(error unknown value "$(REAL_PRECISION)" for REAL_PRECISION)
      endif
    endif
  else
    $(error REAL_PRECISION is not set)
  endif
endif

#### FFTLIB -> _FFTLIB
# allowed input values:   "", "nofft", "FFT_FFTW3", "FFT_MKL", *
# preprocessor variables: "NOFFT", "FFT_FFTW3", "FFT_MKL"
ifneq ($(_FFTLIB),)
  $(error _FFTLIB was passed to make or set via an environment variable)
else
  ifeq ($(FFTLIB),FFT_FFTW3)
    _FFTLIB=FFT_FFTW3
  else ifeq ($(FFTLIB),FFT_MKL)
    _FFTLIB=FFT_MKL
  else ifeq ($(FFTLIB),)
      _FFTLIB=NOFFT
  else	
      _FFTLIB=$(empty)
  endif
endif

#### MPI -> _MPI2
# allowed input values:   "", "mpi", "mpi2", "usempi2", *
# preprocessor variables: "mpi2", "mpi2 mpif90_interface", "nompi2", ""
ifneq ($(_MPI2),)
   $(error _MPI2 was passed to make or set via an environment variable)
else
## Translate mpi (deprecated) into mpi2
## mpi conflicts with the interface module name (issue 87)
  ifeq ($(strip $(MPI)),mpi)
    _MPI2=mpi2
  else
    ifeq ($(strip $(MPI)),mpi2)
      _MPI2=mpi2
    else
## Use the Fortran 90 interface, "use mpi", instead of "include mpif.h"
      ifeq ($(strip $(MPI)),usempi2)
        _MPI2=mpi2 mpif90_interface
      else
        ifneq ($(strip $(MPI)),)
          _MPI2=nompi2
        else
          _MPI2=$(empty)
        endif
      endif
    endif
  endif
endif

#### IMPLICIT -> _IMPLICIT
# allowed input values:   "umfpack", "umfpack32", ""
# preprocessor variables: "umfpack", "umfpack kernel32bit", ""
ifneq ($(_IMPLICIT),)
  $(error _IMPLICIT was passed to make or set via an environment variable)
else
  ifeq ($(IMPLICIT),umfpack32)
    _IMPLICIT=umfpack kernel32bit
  else
    ifneq ($(IMPLICIT),)
      ifeq ($(IMPLICIT),umfpack)
        _IMPLICIT=$(IMPLICIT)
      else
        $(error unknown value "$(IMPLICIT)" for IMPLICIT)
      endif
    endif
  endif
endif

#### PERF -> _PERF
# allowed input values:   * (anything is allowed)
# preprocessor variables: * (matches input)
ifneq ($(_PERF),)
  $(error _PERF was passed to make or set via an environment variable)
else
  _PERF = $(PERF)
endif

#### SLEPC -> _SLEPC
# allowed input values:   "HAVE_SLEPC", *
# preprocessor variables: "HAVE_SLEPC", ""
ifneq ($(_SLEPC),)
  $(error _SLEPC was passed to make or set via an environment variable)
else
  ifeq ($(SLEPC),HAVE_SLEPC)
    _SLEPC=$(SLEPC)
  endif
endif

#### IO_LIB -> _IO_LIB
# allowed input values:   "HAVE_HDF5", *
# preprocessor variables: "HAVE_HDF5", ""
ifneq ($(_IO_LIB),)
  $(error _IO_LIB was passed to make or set via an environment variable)
else
  ifeq ($(IO_LIB),HAVE_HDF5)
    _IO_LIB=$(IO_LIB)
  endif
endif

#### LIBRSB -> _LIBRSB
# allowed input values:   "HAVE_LIBRSB", *
# preprocessor variables: "HAVE_LIBRSB", ""
ifneq ($(_LIBRSB),)
  $(error _LIBRSB was passed to make or set via an environment variable)
else
  ifeq ($(LIBRSB),HAVE_LIBRSB)
    _LIBRSB=$(LIBRSB)
  endif
endif

#### LIBMKL -> _LIBMKL
# allowed input values:   "HAVE_MKL", *
# preprocessor variables: "HAVE_MKL", ""
ifneq ($(_LIBMKL),)
  $(error _LIBMKL was passed to make or set via an environment variable)
else
  ifeq ($(LIBMKL),HAVE_MKL)
    _LIBMKL=$(LIBMKL)
  endif
endif

#### MISC -> _MISC
# allowed input values:   * (anything is allowed)
# preprocessor variables: * (matches input)
ifneq ($(_MISC),)
  $(error _MISC was passed to make or set via an environment variable)
else
  _MISC = $(MISC)
endif

DEFINED_NAMES := $(_REAL_PRECISION) $(_FFTLIB) $(_MPI2) $(_IMPLICIT) $(_PERF) $(_SLEPC) $(_IO_LIB) $(_MISC) $(_LIBRSB) $(_LIBMKL)

C_DEFINED_NAMES := $(_IMPLICIT)

## Usually, preprocessing can be done via -DNAME1=val1 -DNAME2=val2 etc. Set
## preprocessing to work like that here if nothing was specified in the config
## file.
PREPROC_FLAG ?=
PREPROC_SEP ?= $(space)
PREPROC_PREFIX ?=-D

## full preprocessing options - nothing should need changing here.
ifneq ($(DEFINED_NAMES),)
  PREPROCS = $(PREPROC_FLAG)$(subst $(space),$(PREPROC_SEP)$(PREPROC_PREFIX),$(space)$(strip $(DEFINED_NAMES)))
endif

## C preprocessing; like PREPROC_* above.
C_PREPROC_FLAG   ?= $(PREPROC_FLAGS)
C_PREPROC_SEP    ?= $(PREPROC_SEP)
C_PREPROC_PREFIX ?= $(PREPROC_PREFIX)

## GKW C_PREPROCS
ifneq ($(C_DEFINED_NAMES),)
  C_PREPROCS = $(C_PREPROC_FLAG)$(subst $(space),$(C_PREPROC_SEP)$(C_PREPROC_PREFIX),$(space)$(strip $(C_DEFINED_NAMES)))
endif

##############################################################################
### Compiler flag manipulation (see $(TEMPLATE) as defined in the top level  #
### makefile for variable details).                                          #
##############################################################################

ifeq ($(WARN),on)
  ifneq ($(FFLAGS_WARN),)
    FFLAGS += $(FFLAGS_WARN)
  endif
endif

ifneq ($(FFLAGS_INC),)
  FFLAGS += $(FFLAGS_INC)
endif

ifeq ($(DEBUG),on)
  ifneq ($(FFLAGS_DEBUG),)
    FFLAGS += $(FFLAGS_DEBUG) -DDEBUG
  endif
endif

ifeq ($(OPTFLAGS),on)
  ifneq ($(FFLAGS_OPT),)
    FFLAGS += $(FFLAGS_OPT)
  endif
endif

ifeq ($(REAL_PRECISION),real_precision_double)
  ifneq ($(FFLAGS_DOUBLE),)
    FFLAGS += $(FFLAGS_DOUBLE)
  endif
endif

ifeq ($(SMP),OPENMP)
  ifneq ($(FFLAGS_OMP),)
    FFLAGS  += $(FFLAGS_OMP)
  endif
endif

ifneq ($(FFLAGS_OTHER),)
  FFLAGS += $(FFLAGS_OTHER)
endif

## linker LD: same as FC if not specified in config.
ifeq ($(LD),)
  LD = $(FC)
endif

## FFLAGS is added to LD, unless FFLAGS_LD is set in the config.
FFLAGS_LD ?= $(FFLAGS)
LD += $(FFLAGS_LD)

## Set some sensible LDFLAGS for FFTW if nothing is defined.
ifeq ($(FFTLIB),FFT_FFTW3)
  LDFLAGS ?= -lfftw3 -lfftw3f
endif

ifeq ($(SMP),OPENMP)
  LDFLAGS += $(FFLAGS_OMP)
endif

ifeq ($(findstring umfpack,$(IMPLICIT)),umfpack)
  LDFLAGS += $(LD_UMFPACK) $(LD_BLAS)
endif

#APS: workaround for existing config files -- please remove when fixed
LDFLAGS := $(subst lumfpack_gkw32_NB,lumfpack,$(LDFLAGS))
LDFLAGS := $(subst lumfpack_gkw64_NB,lumfpack,$(LDFLAGS))


##############################################################################
### GKW variables (mainly for inclusion in the source code)                  #
##############################################################################

## How to obtain the version

VERSION = $(DEFAULT_GKW_VERSION)
# Number the executable using git describe, when available
# This will work best if tags are added often
VERSION = $(shell git describe --tags --always --dirty 2> /dev/null || echo $(DEFAULT_GKW_VERSION) )

## Fortran compiler variable
ifeq ($(FC),)
  GKW_FC = UNKNOWN
else
  GKW_FC = $(strip $(shell ( type $(FC) || echo "FC_NOT_FOUND" ) | awk '{print $$NF}' ))
endif

## the executable name
EXEC       = $(EXEC_PREFIX)_$(TNAME)_$(VERSION)$(EXEC_SUFFIX)

## input checking executable
CHECKEXEC = $(INPUT_CHECK)_$(TNAME)_$(VERSION)$(EXEC_SUFFIX)

## UMFPACK: library (pre-built or compiled with gkw), BLAS, and the wrapper
## routine.
ifeq ($(findstring umfpack,$(IMPLICIT)),umfpack)
  ifeq ($(UMFPACK_LIBRARY),pre-built)
# no target library if pre-built
  else
    UMF_LIBS = $(UMFPACK_INTERNAL)
    ifneq ($(LD_BLAS),)
      _BLAS=$(LD_BLAS)
    else
      _BLAS=
    endif
  endif
  UMF_INTERFACE=umfpack_wrapper.o
endif

## MKL, for sparse BLAS: it seems that one cannot rely on finding a .mod file
ifeq ($(findstring HAVE_MKL,$(LIBMKL)),HAVE_MKL)
  MKL_SPBLAS_MODULE=mkl_spblas.o
endif

export

##############################################################################
### Targets: 'program' is what creates the gkw executable, provided every    #
### object in '$(OBJLIST)' is built. The version number is used in the       #
### executable name and the included header file 'gkw_info.h'.               #
### The file 'deps.mk' contains a list of object dependencies and should be  #
### automatically updated when any source files are modified.                #
##############################################################################

## Include the list of object files.
include $(SRCDIR)/objfiles.mk

deps: deps.mk

all: program

pre_build: check_config

libs: $(UMF_LIBS) $(MKL_SPBLAS_MODULE)

program: pre_build $(OBJLIST)
	$(LD) -o $(EXEC) $(OBJLIST) $(LDFLAGS)

input_check : pre_build $(CHECKOBJLIST)
	$(LD) -o $(CHECKEXEC) $(CHECKOBJLIST) $(LDFLAGS)

$(UMFPACK_INTERNAL):
	cd $(LIBDIR)/UMFPACK && $(MAKE) umfpack CC='$(CC)' CFLAGS='$(CFLAGS)' CONFIG=''

mkl_spblas.o: $(MKL_SPBLAS_MODULE_SRC)
	echo $(MKL_SPBLAS_MODULE_SRC)
	$(FC) -I$(MKL_HOME)/include $(MKL_HOME)/include/mkl_spblas.f90 -c

## REQUIREDMODULES checking is not ready to use yet
## (and may never be enabled since it can be too restrictive)
## The variable is currently used for error output information only.
# Find missing required modules of the Modules software environment
# management package.
ifeq ($(LOADEDMODULES),)
  _modules_provided := _none_
else
  ifneq ($(REQUIREDMODULES),)
    _modules_provided := $(subst $(colon),$(space),$(LOADEDMODULES))
    ifneq ($(REQUIREDMODULES_SLASH),)
      _modules_provided := $(subst /,$(space)_,$(_modules_provided))
    endif
  else
    _modules_provided :=
  endif
endif
##_modules_missing := $(filter-out $(_modules_provided),$(REQUIREDMODULES))
_modules_missing :=

# Check the Fortran compiler
ifeq ($(FC),)
  _have-Fortran-compiler := not-set
else
  # the word function allows use of FC variables such as "wrapper -options path"
  _have-Fortran-compiler := $(shell command -v $(word 1, $(FC)) > /dev/null || echo "no"  )
endif

## Target which is run before compiling to performs checks.
check_config:
# Some require modules to be loaded
ifneq ($(_modules_missing),)
	$(warning WARNING: module(s) "$(REQUIREDMODULES)" are required.)
  ifeq ($(_modules_provided),_none_)
	  $(warning WARNING: no modules appear to be loaded.)
  else
	  $(warning WARNING: module(s) "$(_modules_missing)" are missing)
  endif
	$(config-error)
endif
# Fortran compiler
ifeq ($(_have-Fortran-compiler),not-set)
	$(warning ERROR: no Fortran compiler FC is set)
	$(config-error)
endif
ifeq ($(_have-Fortran-compiler),no)
	$(warning ERROR: Fortran compiler '$(word 1, $(FC))' NOT FOUND)
	$(config-error)
endif
# usually compiler flags are required to promote double precision
ifeq ($(REAL_PRECISION),real_precision_double)
  ifeq ($(FFLAGS_DOUBLE),)
	  $(warning ERROR: 'FFLAGS_DOUBLE' should be set for REAL_PRECISION=real_precision_double)
	  $(config-error)
  endif
endif
# warn if OMP compiler flags are not set when using OpenMP
ifeq ($(SMP),OPENMP)
  ifeq ($(FFLAGS_OMP),)
	  $(warning WARNING: 'FFLAGS_OMP' should be set for SMP=OPENMP (usually))
	  $(config-error)
  endif
endif
# Big warning if the custom file is not customised
ifneq ($(CONFIG_WARNING),)
	@echo "##### PLEASE EDIT YOUR CONFIG FILE !!! #####"
	@echo "##### PLEASE EDIT YOUR CONFIG FILE !!! #####"
	@echo "##### PLEASE EDIT YOUR CONFIG FILE !!! #####"
	@echo "##### PLEASE EDIT YOUR CONFIG FILE !!! #####"
	@echo "##### PLEASE EDIT YOUR CONFIG FILE !!! #####"
	@echo "      (using $(CONFIGFILE))"
	@echo "##### PLEASE EDIT YOUR CONFIG FILE !!! #####"
	@echo "##### PLEASE EDIT YOUR CONFIG FILE !!! #####"
endif
# (config error message)
define config-error
  $(warning >>> Error occurred when using the file(s))
  $(warning >>>     $(CONFIGFILE))
  $(warning >>>     $(COMPILERFILE))
  $(warning >>> (which may be due to something in that/those file(s)))
  $(warning >>>     It is suggested that you should have these (or similar) modules loaded:)
  $(warning >>>     $(REQUIREDMODULES))
  $(error Cannot build GKW - aborting)
endef

## GKW included information -- only update the file if the contents have
## changed. The create_fortran_file_info() routine can be extended as
## desired to add more information into the source code.
## To ensure the check is always run, a dependency on $(EXEC).info is
## included, which ensures that the executable name and gkw_info.h 
## always agree even after merges and partial rebuilds.
MPIRUNCMD ?= NONE
gkw_info.h: $(EXEC).info $(SRCDIR)/*.*90 $(SRCDIR)/*.mk
	@create_fortran_info_file () \
	{ \
	   echo "character (len=64), parameter, public :: GKW_REV = '$(VERSION)'" \
	         > $$1 ; \
	   echo "character (len=64), parameter, public :: GKW_EXE = '$(EXEC)'"  \
	         >> $$1 ;  \
	   echo "character (len=128), parameter, public :: GKW_FC = '$(GKW_FC)'" \
	         >> $$1 ; \
           echo "character (len=32), parameter, public :: MPIRUNCMD = '$(MPIRUNCMD)'" \
	         >> $$1 ; \
	} ; \
	if [ ! -f "gkw_info.h" ] ; then \
	   create_fortran_info_file "gkw_info.h" ; \
	else \
	   create_fortran_info_file "gkw_info.h.tmp" ; \
	   cmp -s "gkw_info.h" "gkw_info.h.tmp" ; \
	   if [ $$? -ne 0 ] ; then \
	      mv gkw_info.h.tmp gkw_info.h ; \
	   else \
	      rm gkw_info.h.tmp ; \
           fi ; \
	fi

## Target to create an information file. This file might be used for obtaining
## the preprocessor options etc. to recompile the code from a source tarball
## obtained from the executable itself. This would require several makefile
## variables to be output here, into what could be a new makefile. This file
## is related to the gkw_info.h file, so should perhaps be generated in
## conjunction with it. N.B. presently this is a phony target and always run.
%.info: check_config
	@echo "# some makefile variables produced building the code" > $(EXEC).info
	@echo "PREPROCS = $(PREPROCS)" >> $(EXEC).info
	@echo "FFLAGS = $(FFLAGS)" >> $(EXEC).info
	@echo "LDFLAGS = $(LDFLAGS)" >> $(EXEC).info
	@echo "LOADEDMODULES = $(LOADEDMODULES)" >> $(EXEC).info
	@echo "FC = $(FC)" >> $(EXEC).info
	@echo "CONFIG = $(CONFIG)" >> $(EXEC).info
	@echo "EXEC = $(EXEC)" >> $(EXEC).info


## Target to update the module object dependencies file.
deps.mk: $(SRCDIR)/*.*90 $(SRCDIR)/objfiles.mk $(SRCDIR)/gkw.mk $(SRCDIR)/scripts/mkdeps
	@echo "regenerating dependencies"
	@$(SRCDIR)/scripts/mkdeps $(SRCDIR)/objfiles.mk > $(SRCDIR)/deps.mk


##############################################################################
### Compiling rules                                                          #
##############################################################################

%.o: %.F90
	$(FC) $(FFLAGS) $(PREPROCS) -c -o $@ $<

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

%.o: %.c
	$(CC) $(CFLAGS) $(C_PREPROCS) -c -o $@ $<

##############################################################################
### Dependencies                                                             #
##############################################################################

## basic dependencies which should rarely or never change
$(OBJLIST): $(MKFILE) $(CONFIGFILE) $(COMPILERFILE) $(UMF_LIBS)
$(CHECKOBJLIST): $(MKFILE) $(CONFIGFILE) $(COMPILERFILE) $(UMF_LIBS)

## Dependencies for the objects generated from Fortran code; these are updated
## as necessary.
-include $(SRCDIR)/deps.mk
