## --- IFERC helios CSC with eigenvalue solver (eiv) ---
# Inherit defaults
include ${CONFIGDIR}/default.mk

REQUIREDMODULES = intel bullxmpi fftw/3.3 petsc/3.5.3/complex slepc/3.5.3/complex 

#Add settings for eiv libraries petsc and slepc
SLEPC = HAVE_SLEPC
LDFLAGS+= -L${PETSC_DIR}/lib -L${SLEPC_DIR}/lib -lslepc -lpetsc -mkl #-llapack -lblas -lX11
FFLAGS_INC+= -I/${PETSC_DIR}/include -L${SLEPC_DIR}/include
