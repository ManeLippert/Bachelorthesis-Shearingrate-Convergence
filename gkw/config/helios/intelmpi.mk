# inherits settings from default.mk
include ${CONFIGDIR}/default.mk

# to make intelmpi this your personal default, do
# ln -s intelmpi.mk ${USER}.mk 
# and add export USE_IMPI=1 to your login script for gkwnlin
# and load module for intelmpi instead of bullxmpi
REQUIREDMODULES= intel intelmpi fftw/3.3
FC = mpiifort

# Needed for intelmpi I/O, or only for changing striping settings ?
LDFLAGS+=-lirc

