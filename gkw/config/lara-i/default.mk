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
### This file contains a comprehensive list of all makefile variables used by
### GKW.  Simpler examples for specific compilers are in the templates folder.
### Please remove the unnecessary parts of the file before checking into the
### repository. Note that some existing checked in config files may be already
### quite close to what you need and are much more concise than this file.
##############################################################################

## The line below results in a warning, which might be because this file was
## copied directly from the template and not edited. Remove it!
CONFIG_WARNING=yes

## It is possible to split the compiler specific parts of the build
## configuration into a separate file. This is useful when there are multiple
## compilers users may wish to use on a given machine. If COMPILER is set,
## make expects to use the file config/$(HOSTNAME)/compiler/$(COMPILER).mk for
## this purpose. Here $(HOSTNAME) corresponds to a makefile variable which is
## usually derived from the environment variable HOSTNAME. This variable has
## no default value set.
#COMPILER=gnu

## Turn on the optimisation specific Fortran compiler flags which are defined
## in FFLAGS_OPT (see elsewhere in this file). Likely to slow down code
## compilation if enabled, but compilation speed is rarely an issue. One may
## wish to set this to 'off' when DEBUG (see below) is set to 'on' (although
## not necessarily). The default value is 'on'.
#OPTFLAGS=on

## Similar role to OPTFLAGS; turns on flags defined in FFLAGS_DEBUG, which are
## for debugging and code checking purposes. The default value is 'off'. At
## present, this switch does not involve any code preprocessing.
#DEBUG=off

## Similar role to DEBUG; turns on any flags defined in FFLAGS_WARN, which
## are used to produce compile time warnings. The default value is 'on'.
#WARN=on

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
#FFLAGS_WARN=-W -Wall

## The Fortran compiler program (default value is 'gfortran'). When using MPI,
## using a compiler wrapper such as mpif90 can reduce the number or complexity
## of other variables set in this file.
#FC=mpif90

## Compiler debugging flags used if DEBUG=on (default is '-g'). For the
## intel compiler, you might want to put '-g -C -traceback', for example.
#FFLAGS_DEBUG=-g -fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow

## Compiler optimisation flags used if OPTFLAGS=on (default is '-O2')
#FFLAGS_OPT=-O3

## Compiler flag used when compiling in double precision (default is for
## gfortran '-fdefault-real-8 -fdefault-double-8').
## These flags are only used with the corresponding pre-processing setting
## For intel compiler, use -r8
## For IBM compiler, use -qautodbl=dbl4
#FFLAGS_DOUBLE=-r8

## Compiler include flags (default is '-I.')
## If not using a compiler wrapper (such as h5fc or mpif90), appropriate
## FFLAGS_INC may need to be set for including header files required by MPI
## or the FFT library, if not automatically in the include path. The header
## file for FFTW3 is usually installed in /usr/include on GNU/Linux. Note
## that it is also necessary to have the local compilation directory (.) in
## this path to include the compile time information into the code. Most
## compilers use the '-I' flag for include paths to search. Setting FFLAGS_INC
## as in the commented out line below will be sufficient in many cases,
## although may not be necessary.  For some compilers it is necessary to
## explicitly include the present working directory ./
#FFLAGS_INC=-I. -I/usr/include

## OpenMP specific compiler flags; these will only be used if OMP is set in
## preprocessing options below. (-openmp for intel, -fopenmp for gnu)
#FFLAGS_OMP=

## Miscellaneous Fortran compiler flags. One could put, for example,
## architecture-specific flags in here.
## To force a sequential writing of the FDS-file use -DFORCE_IO_ATOMICITY,
## this mitigates an mpi-deadlock which appeared on bt-machines, but could
## slow down the writing process, see issue #251
#FFLAGS_OTHER=-march=native 

## Linker (usually the Fortran compiler). This defaults to FC if not set.
#LD=$(FC)

## By default, all specified FFLAGS_* are used in linking; this can be turned
## off by specifying/uncommenting the FFLAGS_LD line below.
#FFLAGS_LD=

## Flag/variable to enable BLAS in the internally built UMFPACK library (only
## used for implicit and nonspectral runs). If LD_BLAS is non-empty, UMFPACK
## will be built with calls to BLAS routines and the contents of LD_BLAS will
## be added to linking flags.
#LD_BLAS=-lblas

## Variables and flags to link with a pre-built version of UMFPACK (such as that
## usually installed with the suitesparse package). This is only useful if
## implicit or nonspectral runs are required (i.e. 'IMPLICIT=umfpack' is set).
## However, the version of UMFPACK that is provided with GKW is almost always
## sufficient, and the option here is just an alternative for those who can
## not (or do not wish to) build that version. Set 'UMFPACK_LIBRARY=pre-built'
## to use this option. Otherwise, if nothing is set here the GKW provided
## version will be compiled (alternatively, on can set 'UMFPACK_LIBRARY=compile'
## to enforce that). Set LD_UMFPACK to contain the required linking flags for
## the pre-built library. Often, setting only 'LD_UMFPACK=-lumfpack' is
## sufficient, but sometimes extra flags are needed, such as in the examples
## below, where UMFPACK needs to be linked with other parts of suitesparse and
## BLAS (depending on how the library was built).
#UMFPACK_LIBRARY=pre-built
#LD_UMFPACK=-lumfpack
## Suitesparse 4.4.0 with openblas (can usually be switched with -lblas)
#LD_UMFPACK=-lumfpack -lamd -lsuitesparseconfig -lopenblas
## Suitesparse 4.4.0 with CHOLMOD enabled UMFPACK and BLAS
#LD_UMFPACK=-lumfpack -lcholmod -lamd -lsuitesparseconfig -lcolamd -lblas

## Additional linking flags. If no MPI wrapper is used, some flags may need to
## be set here. The OpenMP library could require something. If FFTW3 is used,
## you should probably set at least '-lfftw3 -lfftw3f' to link against double
## and single precision FFTW libraries; these will be set by default if LDFLAGS
## is empty (but this may be incorrect for your system). '-L' is usually used
## to specify directories in which required libraries may be found.
#LDFLAGS=-L/usr/lib -lfftw3 -lfftw3f

### Hint:
### Remember that you can use
###    ldd <executable>
### to see the shared libraries the executable is effectively linked against.
###

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
#PREPROC_FLAG=-WF
#PREPROC_SEP=,
#PREPROC_PREFIX=-D

## If set, the above preprocessor variables may cause the compilation of C 
## code (for umfpack) to fail, if different preprocesor options are required.
## In this case, the C preprocessing can be reset to normal as below
#C_PREPROC_FLAG=
#C_PREPROC_SEP=${space}
#C_PREPROC_PREFIX=-D

## Set the real precision. At present, this value must be set. The default is
## 'real_precision_double', which enables double precision and requires
## the compiler to take all reals in the code as doubles, preferentially via
## the compiler flags in FFLAGS_DOUBLE (which also must be set in this case).
## The other option is 'real_precision_default', which uses the Fortran
## default reals. If using FFTW3, linking must be done with the appropriate
## precision library (see the FFT preprocessing section here).
#REAL_PRECISION=

## Compile MPI functionality into the code. This will require an MPI2
## implementation. Set this to 'mpi2' or 'usempi2', otherwise do not set when
## compiling without an MPI library (or set it to 'nompi'). If not using a
## wrapper (such as h5pfc or mpif90), appropriate FFLAGS_INC and LDFLAGS will
## need to be set.
# Use Fortran 90 interface (is recommended, "use mpi")
#MPI=usempi2
# Use F77 header (used to be more portable, "include mpif.h")
#MPI=mpi2
#MPI=nompi

## Select the FFT library to use. This does not need to be set if FFTs are not
## required (the default). Otherwise, use FFT_FFTW3, then if necessary, add
## something in FFLAGS_INC. With FFTW3 you need to link against different
## library versions, depending on the required precision. It is not a problem
## to link against both, so the commonly used linking flags will be set by
## default if this option is enabled. You may need to set something in LDFLAGS
## if those flags are not correct for your system.
#FFTLIB=FFT_FFTW3

## Select the library used for the implicit and nonspectral schemes.
## This has only been extensively tested on 64 bit in double precision
## But other options should now work.
## For x86_64, select umfpack
## For i686,   select umfpack32
## The default action of setting these flags will be to build the version of
## UMFPACK that is provided with the code. In that case, one requires a C compiler
## which can be set via CC, and corresponding flags via CFLAGS; see the appropriate
## section below. UMFPACK can be built with BLAS support [see LD_BLAS]. Another
## option is to use a pre-built version of UMFPACK [see comments about LD_UMFPACK in
## this file].
#IMPLICIT=umfpack
#IMPLICIT=umfpack32

## Switch, to include the slepc and petsc libraries for the eigenvalue solver.
## If the libraries are not in the standard search path, you must also set
## LDFLAGS and FFLAGS_INC appropriately. You will also need to set -
#SLEPC=HAVE_SLEPC

## Include the HDF5 libraries to write diagnostic output
## and metadata in the HDF5 Format.
#IO_LIB=HAVE_HDF5

## The HDF5 libraries must be in the search path to compile with HDF5
## functionality. This can be achieved by either
##   a) use FC=h5fc to compile with the wrapper provided by the HDF5 library
## or
##   b) setting the LDFLAGS and FFLAGS_INC appropriately. To this end, it
##      may be helpful to copy the correct flags from the
##      output of $(h5fc -show)
##
## In case a), where FC=h5fc, if another compiler wrapper is also required 
## (such as mpif90), you may need to set the h5fc environment variables:
#HDF5_FC=mpif90
#HDF5_FLINKER=mpif90

## Switch, to include the librsb library for fast parallel sparse linear algebra.
## If the libraries are not in the standard search path, you must also set
## LDFLAGS and FFLAGS_INC appropriately.
#LIBRSB = HAVE_LIBRSB

## Enable OpenMP (to use alone or in conjunction with MPI). This may require
## specific Fortran compiler flags to be used. If so, set them in FFLAGS_OMP.
#SMP=OPENMP

## Most can ignore these variables as they only *might* be useful for
## developers testing performance of various parts of the code.
#PERF=perf  # performance timings

## Any additional variables can be passed via MISC. Only really useful for
## code development and testing. It is *strongly* recommended that this should
## *not* be set in any makefile, because no checks are performed on it. It
## appears here for documentation only and if required should be passed directly
## to "make", e.g. make MISC="blas other=somevalue"
#MISC=do_not_set_me

## If compiliaton requires the loading of modules, it can be useful later if
## those needed are listed in the variable below after sucessful compilation.
## This can be obtained from the output of the "module list" command.
## This is for information only, does not currently affect the compilation.
#REQUIREDMODULES=

##############################################################################
### Miscellaneous variables. The items here are listed approximately in
### (descending) order of their perceived usefulness. It is suggested to set
### CC. The others are less likely to be useful.
##############################################################################

## Set the C compiler that can generate binaries for the target machine.
## This is used only to build the provided external libraries if required
## (UMFPACK and AMD, only when implicit and nonspectral runs are desired).
#CC=gcc

## Set the flags to be used with CC. These are just used for optimization
## or debugging of the included external libraries (UMFPACK and AMD) if 
## one needs to override the defaults in the library makefiles.
#CFLAGS=

## The command typically used to run the program with MPI. Examples are
## `mpirun', `mpiexec', `aprun'. At present, this may/will be used to
## facilitate testing of the executable, but very little else. It is set here
## explicitly as there is no obvious straightforward way to gather this
## information at compile time.
#MPIRUNCMD=mpirun

## The first and last parts of the executable name: prefix `gkw' and suffix
## `.x' are the defaults. When compiling _for_ Microsoft Windows, it may be
## useful to set the suffix EXEC_SUFFIX to '.exe'.
#EXEC_PREFIX=gkw
#EXEC_SUFFIX=.exe

## The prefix for the input checking executable name.
#INPUT_CHECK=input_check

## If the code is outside the git version control system, the revision cannot 
## be deduced be the normal method.  It may then be useful to manually set 
## the version used for the executable name.
#DEFAULT_GKW_VERSION=some_version_string

##############################################################################
##############################################################################
