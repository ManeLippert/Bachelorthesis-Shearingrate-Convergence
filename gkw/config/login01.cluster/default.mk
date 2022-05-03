## -- Mesocentre Aix-Marseille University --
# Note: several test cases fail with openmpi (differences in the 5th digit) so use mpich instead
#
# Need to load the following modules before compilation:
# module load userspace/all
# module load gcc
# module load mpich/gcc72/psm2/3.2.1 
# module load fftw3/gcc72/3.3.6-pl2

OPTFLAGS = on
DEBUG = off
PERF = perf

CC = gcc
FC = mpif90  # gfortran
LD = $(FC)

IMPLICIT=umfpack
GKW_LIBS=/home/ycamenen/codes/gkw/libs
MPI = usempi2
#SMP = OPENMP
FFTLIB = FFT_FFTW3

FFLAGS_DOUBLE = -fdefault-real-8 -fdefault-double-8 # gfortran
FFLAGS_OPT = -O3 -ftracer -fomit-frame-pointer -pipe -fweb #gfortran
FFLAGS_OTHER = #-W -Wall -Wunderflow
#FFLAGS_OMP = -fopenmp
FFLAGS_INC = -I. -I$(FFTW3_DIR)/include 
LDFLAGS    = -L$(FFTW3_DIR)/lib/ -lfftw3 -L$(GKW_LIBS)/UMFPACK/Lib/ -lumfpack_gkw64_NB -L$(GKW_LIBS)/AMD/Lib/ -lamd
FFLAGS_DEBUG = -g -fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow #gfortran


##############################################################################
### GKW build configuration template. This file can be used to create a per
### host/user/compiler configuration which will be included by the main GKW
### GNU makefile. All the useful variables that can be set in the makefiles
### should be documented here; please update as necessary.
###
### This is a GNU makefile where lines beginning with '#' are comments.
### Uncomment the described options below as necessary. Where mentioned
### anywhere here that a variable has a default value, this refers to the
### value that was set in the global_defaults.mk file at the time of writing.
### The values in that file will be overridden if specified again here.
###
### Please remove the unnecessary parts of the file before checking into the
### repository. Note that some existing checked in config files may be already
### quite close to what you need and are much more concise than this file.
##############################################################################

## The line below results in a warning, which might be because this file was
## copied directly from the template and not edited. Remove it!
#CONFIG_WARNING = yes

## It is possible to split the compiler specific parts of the build
## configuration into a separate file. This is useful when there are multiple
## compilers users may wish to use on a given machine. If COMPILER is set,
## make expects to use the file config/$(HOSTNAME)/compiler/$(COMPILER).mk for
## this purpose. Here $(HOSTNAME) corresponds to a makefile variable which is
## usually derived from the environment variable HOSTNAME. This variable has
## no default value set.
#COMPILER = gnu

## Turn on the optimisation specific Fortran compiler flags which are defined
## in FFLAGS_OPT (see elsewhere in this file). Likely to slow down code
## compilation if enabled, but compilation speed is rarely an issue. One may
## wish to set this to 'off' when DEBUG (see below) is set to 'on' (although
## not necessarily). The default value is 'on'.
#OPTFLAGS = on

## Similar role to OPTFLAGS; turns on flags defined in FFLAGS_DEBUG, which are
## for debugging and code checking purposes. The default value is 'on'. At
## present, this switch does not involve any code preprocessing.
#DEBUG = on

## If the file is user specific, it might be useful to inherit the defaults
## for the given host. That can be done by uncommenting the following line.
#include $(CONFIGDIR)/default.mk

##############################################################################
### Fortran compiler. You probably need to set FC, FFLAGS_DOUBLE and
### FFLAGS_INC as a minimum. All FFLAGS_* below are collated into a single
### FFLAGS variable at compile time. These flags are also use in linking by
### default, unless the FFLAGS_LD variable is set.
##############################################################################

## Compiler warning flags. If you are modifying the code it is a good idea to
## generate compile-time warning messages.
#FFLAGS_WARN -W -Wall

## The Fortran compiler program (default value is 'gfortran'). When using MPI,
## using a compiler wrapper such as mpif90 can reduce the number or complexity
## of other variables set in this file.
#FC = mpif90

## Compiler debugging flags used if DEBUG = on (default is '-g'). For the
## intel compiler, you might want to put '-g -C -traceback', for example.
#FFLAGS_DEBUG = -g -fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow

## Compiler optimisation flags used if OPTFLAGS = on (default is '-O2')
#FFLAGS_OPT = -O3

## Compiler flag used when compiling in double precision (default is for
## gfortran '-fdefault-real-8 -fdefault-double-8'). 
## These flags are only used with the corresponding pre-processing setting
#FFLAGS_DOUBLE = -r8

## Compiler include flags (default is '-I.')
## If not using a compiler wrapper (such as h5pfc or mpif90), appropriate
## FFLAGS_INC may need to be set for including header files required by MPI
## or the FFT library, if not automatically in the include path. The header
## file for FFTW3 is usually installed in /usr/include on GNU/Linux. Note
## that it is also necessary to have the local compilation directory (.) in
## this path to include the compile time information into the code. Most
## compilers use the '-I' flag for include paths to search. Setting FFLAGS_INC
## as in the commented out line below will be sufficient in many cases,
## although may not be necessary.
#FFLAGS_INC = -I. -I/usr/include

## OpenMP specific compiler flags; these will only be used if OMP is set in
## preprocessing options below.
#FFLAGS_OMP =

## Miscelaneous Fortran compiler flags. One could put, for example,
## architecture-specific flags in here.
#FFLAGS_OTHER = -march=native

## Linker (usually the Fortran compiler). This defaults to FC if not set.
#LD = $(FC)

## By default, all specified FFLAGS_* are used in linking; this can be turned
## off by specifying/uncommenting the FFLAGS_LD line below.
#FFLAGS_LD =

## Additional linking flags. If no MPI wrapper is used, some flags may need to
## be set here. The OpenMP library could require something. If FFTW3 is used,
## you should probably set at least '-lfftw3 -lfftw3f' to link against double
## and single precision FFTW libraries; these will be set by default if LDFLAGS
## is empty (but this may be incorrect for your system). '-L' is usually used
## to specify directories in which required libraries may be found.
#LDFLAGS= -L/usr/lib -lfftw3 -lfftw3f

##############################################################################
### Options that affect Fortran source preprocessing. These variables are
### usually filtered in the makefile to avoid passing unusual variables on
### to the compiler, which may produce less informative error messages (if
### any).
##############################################################################

## With most Fortran compilers, preprocessing can be done by passing
## -D_NAME1=VAL1 -D_NAME2=VAL2 etc. to the compiler.  However, some compilers,
## such as the IBM XL Fortran compiler, require something like
## -WF,-D_NAME1=VAL1,-D_NAME2=VAL2 etc. In that case, the variables below can
## be enabled.
#PREPROC_FLAG   = -WF
#PREPROC_SEP    = ,
#PREPROC_PREFIX = -D

## Set the real precision. At present, this value must be set. The default is
## 'real_precision_double', which enables double precision and requires
## the compiler to take all reals in the code as doubles, preferentially via
## the compiler flags in FFLAGS_DOUBLE (which also must be set in this case).
## The other option is 'real_precision_default', which uses the Fortran
## default reals. If using FFTW3, linking must be done with the appropriate
## precision library (see the FFT preprocessing section here).
#REAL_PRECISION =

## Compile MPI functionality into the code. This will require an MPI
## implementation. Set this to 'mpi' to use and otherwise do not set when
## compiling without an MPI library (or set it to 'nompi'). If not using a
## wrapper (such as h5pfc or mpif90), appropriate FFLAGS_INC and LDFLAGS will
## need to be set.
#MPI = mpi2

## Select the FFT library to use. This does not need to be set if FFTs are not
## required (the default). Otherwise, use FFT_FFTW3, then if necessary, add
## something in FFLAGS_INC. With FFTW3 you need to link against different
## library versions, depending on the required precision. It is not a problem
## to link against both, so the commonly used linking flags will be set by
## default if this option is enabled. You may need to set something in LDFLAGS
## if those flags are not correct for your system.
#FFTLIB = FFT_FFTW3

## Select the external library used for the implicit scheme. The only option
## at present is  'umfpack', which requires the libraries to be built. That
## option is only know to be working with a subset of code options when using
## double precision on x86_64.
#IMPLICIT = umfpack

## Enable OpenMP (to use alone or in conjunction with MPI). This may require
## specific Fortran compiler flags to be used. If so, set them in FFLAGS_OMP.
#SMP = OPENMP

## Most can ignore these variables as they only *might* be useful for
## developers testing performance of various parts of the code.
#PERF = perf  # performance timings

## If compiliaton requires the loading of modules, it can be useful later if 
## those needed are listed in the variable below after sucessful compilation. 
## This can be obtained from the output of the "module list" command.
## This is for information only, is does not currently affect the compilation.
#REQUIREDMODULES = 

##############################################################################
### Miscellaneous variables. The items here are listed approximately in
### (descending) order of their perceived usefulness. It is suggested to set
### CC (and possibly MPIRUNCMD if applicable). The others are less likely to
### be useful.
##############################################################################

## Set the C compiler that can generate executables able to run on the code
## compiling machine (rather than on the target machine). At present, it is
## not necessary to have a C compiler at all; it is only used to compile a
## small C program which can be run to obtain the code revision before
## compiling the Fortran source code. You may wish set something here if after
## building the code your executable name does not appear to contain any
## version information, which may be of some use.
#CC = gcc

## The command typically used to run the program with MPI. Examples are
## `mpirun', `mpiexec', `aprun'. At present, this may/will be used to
## facilitate testing of the executable, but very little else. It is set here
## explicitly as there is no obvious straightforward way to gather this
## information at compile time.
#MPIRUNCMD = mpirun

## The first and last parts of the executable name: prefix `gkw' and suffix
## `.x' are the defaults. When compiling _for_ Microsoft Windows, it may be
## useful to set the suffix EXEC_SUFFIX to '.exe'.
#EXEC_PREFIX = gkw
#EXEC_SUFFIX = .exe

## The prefix for the input checking executable name.
#INPUT_CHECK=input_check

## In the case no C compiler is found (see the CC variable) and the code
## revision can not be deduced be the normal method, it may be useful manually
## set the version used for the executable name.
#DEFAULT_GKW_VERSION= some_version_string

##############################################################################
##############################################################################
