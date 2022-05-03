#############################################################################
###   GKW top-level makefile
###
###
###   With this makefile, GKW will be built using a specific configuration
###   file appropriate for the local machine/user/compilers. These files are
###   found in the config/ directory and are included as necessary or
###   requested into this makefile and the main makefile, src/gkw.mk.
###
###   To build the code, you should not need to edit this file or any file
###   with extension .mk found in the source directory src/ -- it should be
###   possible to add new files to the config/ directory or specify an
###   existing one to use. It is possible to go ahead and use the defaults,
###   but this may not be very useful unless you require only the basic code
###   functionality. A comprehensive example is given below.
###
###   Within the config/ directory, the file "global_defaults.mk" contains
###   variables which are always used, unless overridden by a machine specific
###   file, which would usually be found in a subdirectory corresponding to
###   the machine name. You may wish to edit something in there if, for
###   example, you will always run in single precision.
###
###   The file config/template.mk is an example template file which should
###   document all the variables that might be useful.
###
###   To setup on a new machine, you can try 'make config-file', which will
###   merely copy a template file into the appropriate directory. This file
###   will then likely need to be edited. The name of the file should be
###   "default.mk", or alternatively "${USER}.mk", where ${USER} is your
###   username. This allows for a particular user to override the variables
###   specified in "default.mk". Instead of copying the template file, it
###   may be easier to find and modify an existing file which approximately
###   matches your machine specification.
###
###   The names of the subdirectories of config/ are based on the machine
###   host name (which usually comes from the environment variable HOSTNAME.)
###   This makefile tries to match the directory name with the local value of
###   $(HOSTNAME). If no directory name matches $(HOSTNAME), a directory with
###   at least the first 3 letters in common with $(HOSTNAME) will be
###   considered appropriate. This is useful when 1 file is required for
###   multiple similarly named login nodes.
###
###   To attempt to build with an existing configuration file for another
###   machine, you can try
###      > make HOSTNAME=<some-hostname>
###   which is just really pretending to be that machine. Alternatives are
###      > make CONFIGFILE=<specific-filename>
###   which requires the full filename to use,
###      > make CONFIG=<relative-filename>
###   which selects a file relative to the location of this makefile, e.g.
###   'make CONFIG=config/defaults.mk' would likely be valid. Also
###      > make TARGET=<machine-name>
###   will use config/$(TARGET)/default.mk
###   Additionally, you can do
###      > make USER=some-other-user
###   which will (usually) use config/$(HOSTNAME)/some-other-user.mk
###
###   Because some will want to use different compilers on the same machine,
###   there is the option to put compiler specific settings in a
###   subdirectory 'compiler' of $(HOSTNAME), which is specified via
###   $(COMPILER). This can be set in the primary included makefile
###   $(CONFIGFILE) or when running make e.g. 'make COMPILER=gnu', which will
###   force make to look for the file "compiler/gnu.mk" inside the host
###   specific directory. It is recommended that the compiler directory is
###   used, particularly when multiple users will access that machine.
###
###   If you can not or do not want to use this file to build the code, a much
###   simpler makefile is provided in the $(SRCDIR) (typically src/makefile) -
###   see that file for details.
###
###   As a fairly general example, consider 2 similar machines "vip01" and
###   "vip02", together with 2 users "user1" and "user2", and  2 compilers,
###   "gnu" and "intel". The defaults for both these machines can be set in
###   the file "vip/default.mk". That file can contain the line
###
###     COMPILER=gnu
###
###   which instructs make to include the settings for the gnu compiler
###   contained in "vip/compiler/gnu.mk". Now user2 may wish to use "intel"
###   as the default compiler. In that case, an additional file "vip/user2.mk"
###   can be made and should contain the following 2 lines:
###
###     include $(CONFIGDIR)/default.mk
###     COMPILER=intel
###
###   this means that user2 inherits all of the defaults in the "default.mk"
###   file, but uses instead "vip/compiler/intel.mk" for specific
###   compiler options.
###
###   To build an executable that checks the GKW input file format:
###
###     make checker
###
###   (on hosts which cross-compile for different node hardware, building a
###   checker for the host may require different compiler options from the
###   main executable, e.g. "make clean; make checker MPI=no_mpi FC=gfortran")
###
##############################################################################

##############################################################################
### Set the configuration specific variables via hostname, username and any  #
### variables passed to make. These include the source+build directory,      #
### executable names and makefile names. Further makefiles to be included    #
### depend on some of the values selected here.                              #
##############################################################################

## shell (can be set to something simple for testing purposes)
#SHELL = /usr/bin/dash

## executable filename prefix and suffix
EXEC_PREFIX = gkw
EXEC_SUFFIX = .x

## basic directories
BASEDIR = $(shell pwd)
SRCDIR     = $(BASEDIR)/src
RUNDIR     = $(BASEDIR)/run
MAINOBJDIR = $(BASEDIR)/obj
CONFDIR    = $(BASEDIR)/config
LIBDIR     = $(BASEDIR)/libs

## GKW_HOME is used in some config files
ifeq ($(GKW_HOME),)
  GKW_HOME = $(BASEDIR)
endif

## various makefiles
MKFILE     = $(SRCDIR)/gkw.mk
DEFAULTS   = $(CONFDIR)/global_defaults.mk
TEMPLATE   = $(CONFDIR)/template.mk

## Set HOSTNAME and USER if not provided. LOGNAME is not the same as USER
## but should usually be sufficient.
HOSTNAME ?= $(strip $(shell uname -n || hostname ))
USER     ?= $(strip $(shell [ ! "${LOGNAME}" == "" ] && echo "${LOGNAME}" || \
              logname 2>&1 | grep -v "no login name" || echo "default"))

## Create a list of directory names to be tried for a configuration file.
## e.g. for HOSTNAME=localhost.localdomain, the list would be
## localhost.localdomain, localhost, localhos, localho, ... , loc.
## (Is there an easier way to do this using basic POSIX functionality?)
HLIST = $(HOSTNAME) $(strip $(shell HN=$$(echo $(HOSTNAME) | \
	  sed -e 's/\..*$$//'); until [ $${\#HN} -lt  3 ] ; do echo $$HN ; \
          LHN=$${\#HN} ; HN=$$(echo $$HN | awk -v var=$$LHN '{ string=substr($$1, 1, var - 1); print string; }') ; done ))

## Specify the configuration file relative to the present directory via CONFIG
## or TARGET; if not set here it will be set below.
ifneq ($(CONFIG),)
  CONFIGFILE = $(BASEDIR)/$(CONFIG)
endif
ifneq ($(TARGET),)
  CONFIGFILE = $(CONFDIR)/$(TARGET)/default.mk
else
  TARGET = default
endif

## Try (if not set already) to find a config file using the list of
## directories in HLIST; prefer $(USER).mk to default.mk and default to
## the file $(DEFAULTS) in the top level if one cannot be found.
CONFIGFILE ?= $(strip $(shell for fileprefix in $(USER) default ; do \
             for name in $(HLIST); do \
                 if [ -e $(CONFDIR)/$${name}/$${fileprefix}.mk ]; then \
                     echo "$(CONFDIR)/$${name}/$${fileprefix}.mk"; \
                     found="y" ; \
                     break 2 ; \
                 fi ; \
             done ; done ; \
             if [ ! "$${found}" = "y" ]; then echo "$(DEFAULTS)" ; fi))

## block the use of defaults unless USE_DEFAULTS=true
ifeq ($(CONFIGFILE),$(DEFAULTS))
  ifneq ($(USE_DEFAULTS),true)
    NOCONFIGFILE = true
  else
    CONFIGNAME = default
  endif
endif

## Do we have the primary configuration file? (should have, since $(DEFAULTS) is typically present)
_have_config_file := $(strip $(shell test -e $(CONFIGFILE) && echo "yes" || echo "no"))

## A name derived from the configuration file to compare with TARGET below
CONFIGNAME ?= $(subst $(CONFDIR)/,,$(strip $(shell dirname $(CONFIGFILE))))

# The CONFIGDIR variable is used to point to the directory in which the
# configuration file is present. It is useful in including files within
# makefiles which reside in that directory.
CONFIGDIR ?= $(strip $(shell dirname $(CONFIGFILE)))

## shorted hostname (remove everything after first ".")
HN = $(strip $(shell echo $(HOSTNAME) | sed -e 's/\..*$$//' ))

NOCONFIGFILE ?= false

##############################################################################
### Include the various configuration makefiles: $(DEFAULTS), the            #
### hostname/user specific file and a compiler specific file if set.         #
##############################################################################

include $(DEFAULTS)

ifeq ($(_have_config_file),yes)
  include $(CONFIGFILE)
endif

ifneq ($(COMPILER),)
  COMPILERFILE = $(strip $(shell echo `dirname $(CONFIGFILE)`/compiler/$(COMPILER).mk))
  _have_compiler_file := $(strip $(shell [ -e $(COMPILERFILE) ] && echo "yes" || echo "need"))
else
  _have_compiler_file := no
endif

ifeq ($(_have_compiler_file),yes)
  include $(COMPILERFILE)
endif

##############################################################################
### Based on REAL_PRECISION as set in the configuration files, set the build #
### directory suffix so that single and double precision executables can be  #
### updated independently. Set suffixes for other configuration options.     #
##############################################################################

ifeq ($(REAL_PRECISION),real_precision_double)
  OBJDIR_SUFFIX=-DP
else
  OBJDIR_SUFFIX=-SP
endif

ifeq ($(DEBUG),on)
  DEBUG_SUFFIX=-DEBUG
endif

ifeq ($(OPTFLAGS),off)
  OPT_SUFFIX=-no_opt
endif

## Further granulation of object directories?
#ifneq ($(MPI),mpi)
#  MPI_SUFFIX=-no_MPI
#endif
#
#ifeq ($(SMP),openmp)
#  OMP_SUFFIX=-omp
#endif
#
#ifeq ($(FFTLIB),)
#  FFT_SUFFIX=-no_FFT
#endif

##############################################################################
### Set TNAME, used mainly for the object directory and executable names.    #
##############################################################################

ifneq ($(COMPILER),)
  TNAME_SUFFIX=-$(COMPILER)$(OBJDIR_SUFFIX)
else
  TNAME_SUFFIX=$(OBJDIR_SUFFIX)
endif

ifneq ($(COMPILER),)
  TDIRNAME_SUFFIX=-$(COMPILER)$(OBJDIR_SUFFIX)$(DEBUG_SUFFIX)$(MPI_SUFFIX)$(OMP_SUFFIX)$(OPT_SUFFIX)
else
  TDIRNAME_SUFFIX=$(OBJDIR_SUFFIX)$(DEBUG_SUFFIX)$(MPI_SUFFIX)$(OMP_SUFFIX)$(OPT_SUFFIX)
endif



ifeq ($(TARGET),default)
  ifeq ($(_have_config_file),yes) # call by config name
    TNAME    = $(CONFIGNAME)$(TNAME_SUFFIX)
    TDIRNAME = $(CONFIGNAME)$(TDIRNAME_SUFFIX)
  else # just use the host name
    TNAME    = $(HN)$(TNAME_SUFFIX)
    TDIRNAME = $(HN)$(TDIRNAME_SUFFIX)
  endif
else
  TNAME    = $(TARGET)$(TNAME_SUFFIX)
  TDIRNAME = $(TARGET)$(TDIRNAME_SUFFIX)
endif

##ifeq ($(CONFIGFILE),$(DEFAULTS))
##  TNAME    = $(HN)$(TNAME_SUFFIX)
##  TDIRNAME = $(HN)$(TDIRNAME_SUFFIX)
##else
##  TNAME    = $(TARGET)$(TNAME_SUFFIX)
##  TDIRNAME = $(TARGET)$(TDIRNAME_SUFFIX)
##endif

##############################################################################
### set the object directory, OBJDIR                                         #
##############################################################################

OBJDIR = $(MAINOBJDIR)/$(TDIRNAME)

##############################################################################
### exports to further makefiles                                             #
##############################################################################

export MKFILE CONFIGDIR CONFIGFILE COMPILERFILE DEFAULTS SRCDIR LIBDIR
export EXEC_PREFIX EXEC_SUFFIX TNAME COMPILER INPUT_CHECK GKW_HOME

##############################################################################
### targets (gkw is the default target)                                      #
##############################################################################

all: gkw

## Reports if the requested/required included makefiles are present. By
## default, stop before compiling if the files are missing.
check_target:
ifeq ($(_have_config_file),no)
	$(warning ERROR:)
	$(warning >>> Primary configuration file)
	$(warning >>> $(CONFIGFILE))
	$(warning >>> NOT FOUND.)
	$(warning >>> Try)
	$(warning >>> $(MAKE) $(MAKEFLAGS) config-file)
	$(warning >>> to create a template config file)
	$(warning >>>   *OR* alternatively)
	$(warning >>> Specify the primary configuration file directly via)
	$(warning >>> $(MAKE) $(MAKEFLAGS) CONFIGFILE=<config-file>)
	$(warning >>> (A primary configuration file should be named $(USER).mk)
	$(warning >>>  or default.mk and is found in a subdirectory of)
	$(warning >>> $(CONFDIR).)
	$(error cannot continue)
else
  ifeq ($(NOCONFIGFILE),true)
	  $(warning ERROR:)
	  $(warning >>> NO CONFIGURATION FILE FOUND!)
	  $(warning >>>)
	  $(warning >>> --- PLEASE CREATE A CONFIGURATION ---)
	  $(warning >>> in the file)
	  $(warning >>> $(CONFDIR)/$(HOSTNAME)/default.mk)
	  $(warning >>> or)
	  $(warning >>> $(CONFDIR)/$(HOSTNAME)/$(USER).mk)
	  $(warning >>> You can copy an existing file from one)
	  $(warning >>> of the subdirectories of)
	  $(warning >>> $(CONFDIR))
	  $(warning >>> or you can try)
	  $(warning >>> $(MAKE) $(MAKEFLAGS) config-file)
	  $(warning >>> to create a template config file.)
	  $(warning >>>)
	  $(warning >>> --- ALTERNATIVELY ---, )
	  $(warning >>> specify the primary configuration file directly via)
	  $(warning >>> $(MAKE) $(MAKEFLAGS) CONFIGFILE=<config-file>)
	  $(warning >>>)
	  $(warning >>> --- INSTEAD ---, )
	  $(warning >>> to use the defaults, pass)
	  $(warning >>> USE_DEFAULTS=true)
	  $(warning >>> to make (good luck!))
	  $(error cannot continue)
  endif
endif
ifeq ($(_have_compiler_file),need)
	$(warning ERROR:)
	$(warning >>> Requested compiler configuration file)
	$(warning >>> $(COMPILERFILE))
	$(warning >>> NOT FOUND.)
	$(warning >>> The variable)
	$(warning >>> COMPILER=$(COMPILER))
	$(warning >>> has been set via)
	$(warning >>> $(CONFIGFILE))
	$(warning >>> OR when running make, )
	$(warning >>> OR via the environment variable COMPILER.)
	$(warning >>> Please fix!)
	$(error cannot continue)
endif
ifeq ($(_have_config_file),yes)
	@echo ""
	@echo "Attempting to build with"
	@echo "   $(CONFIGFILE)"
  ifeq ($(_have_compiler_file),yes)
	@echo "and compiler file"
	@echo "   $(COMPILERFILE)"
  endif
	@echo ""
endif


## Generate a template file for a specific machine containing most of the
## useful variables and their descriptions.
config-file:
	@make_template_file() \
	{ \
	    cat $(TEMPLATE) > $$1; \
	    echo ">>> Created a new config file from template." ; \
	    echo ">>> Please now edit the file:"; echo "";  \
	    echo "      $$1," ; echo "" ; \
	    echo ">>> otherwise basic defaults (no FFT or MPI)" ;\
	    echo ">>> can now be used by re-running make..." ; \
	} ; \
	if [ "$(COMPILER)" = "" ]; then \
	  configfile="$(CONFDIR)/$(HOSTNAME)/default.mk"; \
	  test -d $(CONFDIR)/$(HOSTNAME) || mkdir $(CONFDIR)/$(HOSTNAME); \
	else \
	  configfile="$(CONFDIR)/$(HOSTNAME)/compiler/$(COMPILER).mk"; \
	  test -d $(CONFDIR)/$(HOSTNAME) || mkdir $(CONFDIR)/$(HOSTNAME); \
	  test -d $(CONFDIR)/$(HOSTNAME)/compiler || mkdir $(CONFDIR)/$(HOSTNAME)/compiler; \
	fi; \
	if [ ! -f "$${configfile}" ]; then \
	    make_template_file $${configfile} ; \
	else \
	    echo "FILE $${configfile} EXISTS"; \
	    echo "(remove the file to recreate from template)"; \
	fi

## Create various directories for building the code and storing the
## executables.
setup: check_target
	@test -d $(RUNDIR) || mkdir $(RUNDIR)
	@test -d $(MAINOBJDIR) || mkdir $(MAINOBJDIR)
	@test -d $(OBJDIR) || mkdir $(OBJDIR)

## At present, the code is always re-linked when make is run, since the
## executable is moved into the run directory afterwards.
gkw: setup deps
	$(MAKE) -C $(OBJDIR) -f $(MKFILE) all
	cp --preserve=timestamps -f $(OBJDIR)/*$(EXEC_SUFFIX) $(OBJDIR)/*.info $(RUNDIR)

## Build the input checker.
checker: setup deps
	$(MAKE) -C $(OBJDIR) -f $(MKFILE) input_check
	cp --preserve=timestamps -f $(OBJDIR)/*$(EXEC_SUFFIX) $(OBJDIR)/*.info $(RUNDIR)

## Clean the (host,compiler,precision) specific build directory
cleanloc:
	-rm -f $(OBJDIR)/*

## Clean *all* the build directories
clean:
	-rm -fr $(MAINOBJDIR)
	-rm -fr $(SRCDIR)/deps.mk

## Make libraries if required.
libs:
	$(MAKE) -C $(SRCDIR) -f $(MKFILE) libs

## Clean *all* the lib directories
cleanlibs:
	cd $(LIBDIR)/AMD && $(MAKE) purge
	cd $(LIBDIR)/UMFPACK && $(MAKE) purge

## Update the Fortran dependencies.
deps:
	$(MAKE) -C $(SRCDIR) -f $(MKFILE) deps

## build (some of) the documentation, requires a latex installation
docs:
	$(MAKE) -C doc all && echo "Main documentation can now be found in doc/doc.pdf"

## Print some of this file to screen in case you can't just look at it!
help:
	@echo
	@head -n 99 GNUmakefile

tests:
	./scripts/gkw_run_tests

# Do not try to build files with following names, execute target instead
.PHONY: tests help docs deps clean cleanloc config-file libs cleanlibs

##############################################################################
