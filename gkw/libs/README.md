
# Licensing

This directory contains third party libraries that are used by GKW.
These libraries have their own licences that can be different from the 
GKW license. For more information see the README or LICENCE files in 
the corresponding directories. 

# Purpose

UMFPACK is a library used by the implicit and nonspectral schemes.

AMD is a library needed by the UMFPACK library.

The folder UFconfig contains configuration files needed to build UMFPACK.

# Installation

These libraries are compiled automatically during the gkw build process.

    gkwmake -j

They can also be built manually using

    gkwmake libs

To clean them up and delete old compiled files you may want to use the target

    gkwmake cleanlibs

All required makefile settings for UMFpack 
should now be passed through the GKW config files. 

To build the umfpack libraries and link them to GKW, it should in most 
cases be sufficient to set, in the GKW config file:

    IMPLICIT = umfpack  (x86_64)

or

    #IMPLICIT = umfpack32 (i686)

On machines with 32bit kernels (i686),   use the 32bit version of UMFpack.
On machines with 64bit kernels (x86_64), use the 64bit version of UMFpack.

In some cases it may also be necessary to set:

    XCC = gcc   # C compiler for building umfpack

External UMFPACK libraries can also be used, for these and other linking
options, please refer to the comprehensive GKW makefile documentation in 
../config/template.mk.

A note on the deprecated UMFPACK build process:

The scripts compile_umf and compile_umf32 are deprecated, the GKW
make process should now be used.  The UMFPACK Fortran demo and wrapper 
no longer needs to be built, as GKW now has its own C UMFPACK wrapper,
in gkw/src.  Old config files on some machines may still need updating 
to the new process, so the old compile scripts remain for now.

See issue 117 at
https://bitbucket.org/gkw/gkw/issues/117
for more information.
