## --- IFERC helios CSC with serial hdf5 ---
# Inherit defaults
include ${CONFIGDIR}/default.mk

IO_LIB = HAVE_HDF5
REQUIREDMODULES = intel bullxmpi fftw/3.3 hdf5 
HDF5_FLINKER=mpif90
HDF5_FC=mpif90
FC=h5fc
#OMPI_FC=h5fc
#FC=mpif90

