##
##	PDP1 Makefile 
##
FILE_EXT = 
##
# note that Xgrafix must be double precision
XGPATH = /usr/local
# sometimes only version-specific libtk/tcl8.x are in libdirs, if so uncomment:
# TKV=8.4
TKV=8.6
## for macs, use MacPorts Tcl/Tk:
# TKDIR=-L/opt/local/lib
##	Gas type for the mcc package.  Uncomment only one.
##
##
##
EXEC = ../xpdp1
##
CC = gcc
##
##      C compiler.  Normally gcc.
##
CFLAGS= -O3 -I$(XGPATH)/include -DPDP1PATH=\"$(HOME)/ptsg_codes/xpdp1/src\" -Wno-unused-result
#
## CFLAGS= -g -DDEBUG -I\usr\local\include -I$(XGPATH)/include -DPDP1PATH=\"$(HOME)/xpdp1/src\"
##
##      Flags used for compiling.  -O2 for optimization, which is optional.
##	If the X11 include files are not located in /usr/include/X11, you
##	must specify the directory as a flag prefixed by -I.  For instance, 
##	if the include files for X11 are located in /usr/local/include, the 
##	flag -I/usr/local/include must also be passed.
##
LIBS  = -L$(XGPATH)/lib -L/usr/local/lib  -L/usr/lib/x86_64-linux-gnu -lXGC $(TKDIR) -ltk$(TKV) -ltcl$(TKV) -lXpm -lX11 -lm -ldl 
##
## 	Libraries and their directories used in loading.  Normal is -lm and -lX11
##	for the Math and X11 libraries, which is NOT optional.  On Unicos, -lnet
##	is also used.  If the X11 libraries are not located in /usr/lib/X11, you
##	must also specify the directory as a flag prefixed by -L.  For instance,
##	if the include files for X11 are located in /usr/local/lib, the flag
##	-L/usr/local/lib must also be passed.  On Sun, use -Bstatic to avoid
##	problems with differing versions of OS.  Man cc for more information on flags.
##
##
PDP1OBJ= fft.o field.o move.o gather.o pdp1.o start.o load.o prest.o \
	 padjus.o initwin.o maxwellv.o  mccdiaginit.o xsect.o\
	 argonmcc.o heliummcc.o neonmcc.o oxygenmcc.o mcc.o

all:	$(PDP1OBJ) $(EXEC)

.c.o:	pdp1.h xsect.h
	$(CC) -c $(CFLAGS) $*.c

$(EXEC):	$(PDP1OBJ)
		$(CC) -o $(EXEC) $(PDP1OBJ) $(LIBS)

clean:
	@rm *.o 
