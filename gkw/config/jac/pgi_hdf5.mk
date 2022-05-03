## --- JET JAC, pgi with serial hdf5 ---

# Inherit defaults
include ${CONFIGDIR}/pgi.mk

IO_LIB = HAVE_HDF5
REQUIREDMODULES = pgi openmpi hdf5/1.8.5
FC=mpif90

# note that the default module hdf5/1.8.9
# does not appear to have the fortran .mod files
# override this using older version
HDF5_HOME=/usr/local/depot/hdf5-1.8.5/

# h5fc wrpapper can't wrap mpif90, so don't use
LDFLAGS+= -L/${HDF5_HOME}/lib -lhdf5 -lhdf5_fortran -lhdf5_hl -lhdf5hl_fortran -lz -lm
FFLAGS_INC+= -I/${HDF5_HOME}/include

