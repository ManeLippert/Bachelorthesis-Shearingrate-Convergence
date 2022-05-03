PROG          = chease
PROG_HDF5     = chease_hdf5
PROG_ITM      = chease_itm
PROG_ITM_HDF5 = chease_itm_hdf5
PROG_KEPLER   = chease_kepler
PROG_interpos   = chease_interpos

# SRCS_0 chease subroutine in chease.f90
SRCS_0 =	a_chease.f90 chease.f90 acopy.f90 aldlt.f90 apcoef.f90 apcoef2.f90 \
	assign_code_parameters.f90 atcoef.f90 auxval.f90 \
	away.f90 ballit.f90 baloon.f90 basis1.f90 basis2.f90 basis3.f90 \
	basis4.f90 blines.f90 bltest.f90 bndspl.f90 bound.f90 bsexpeq.f90 \
	bsfunc.f90 bstnzpro.f90 ccopy.f90 center.f90 check.f90 \
	chipsi.f90 chipsimetrics.f90 cint.f90 conver.f90 copyap.f90 copyapp.f90 \
	copyat.f90 cotrol.f90 cubrt.f90 curent.f90 cvzero.f90 data.f90 \
	direct.f90 drhodp.f90 dwy.f90 energy.f90 eqdim.f90 erdata.f90 \
	errorch.f90 evlate.f90 four1.f90 fourfft.f90 fourier.f90 g_0.f90 \
	g_1.f90 g_2.f90 g_3.f90 gauss.f90 gchi.f90 gdataext.f90 genout.f90 \
	gijlin.f90 gloadd.f90 globals_init.f90 gloqua.f90 \
	guess.f90 iarray.f90 identa.f90 identb.f90 indexx.f90 initia.f90 \
	iodisk.f90 isamin.f90 ismax.f90 ismin.f90 isofind.f90 \
	isofun.f90 isrchfge.f90 issum.f90 itipr.f90 ivar.f90 jnovaw.f90 \
	labrun.f90 limita.f90 limitb.f90 ltxw.f90 lyv.f90 magaxe.f90 \
	mappin.f90 matrix.f90 mesage.f90 mesh.f90 msplcy.f90 \
	mspline.f90 nerat.f90 nonlin.f90 norept.f90 ntridg.f90 oarray.f90 \
	oldeq.f90 oldnew.f90 outgload.f90 outmetric.f90 outmksa.f90 outnvw.f90 outpen.f90 \
	output.f90 outxt.f90 packme.f90 packmep.f90 page.f90 polyfun.f90 polynm.f90 \
	ppbstr.f90 pprime.f90 pprm.f90 ppspln.f90 ppspln2.f90 premap.f90 \
	preset.f90 prfunc.f90 priqqu.f90 prnorm.f90 profile.f90 psibox.f90 \
	psicel.f90 psvol.f90 qplacs.f90 rarray.f90 realft.f90 reseti.f90 \
	resetr.f90 resppr.f90 rmrad.f90 rscale.f90 runtim.f90 rvar.f90 rvar2.f90 \
	rzbound.f90 saxpy.f90 scopyr.f90 setupa.f90 \
	setupb.f90 shave.f90 smooth.f90 solovev.f90 solvit.f90 sort3.f90 \
	splcy.f90 splcyp.f90 splifft.f90 spline.f90  ssum.f90 \
	stchps.f90 stepon.f90 subsz.f90 surface.f90 surfadd.f90 surfrz.f90 \
	tcase.f90 test.f90 tetare.f90 tpsi.f90 tricyc.f90 tricycm.f90 \
	tridagm.f90 tshift.f90 vacufft.f90 vacuum.f90 vlion.f90 vzero.f90 \
	whtext.f90 witext.f90 wrtext.f90 wrtplot.f90 wrtmat.f90 wrtbin.f90 \
	hamada.f90 neoart.f90 outgyro.f90 bscoeff.f90 outelit.f90 outastro.f90 ogyropsi.f90

SRCS_0_in_interpos = sscal.f90 scopy.f90 sdot.f90 

# hdf5 interface for outputs
SRCS_0_NOHDF5 = write_ogyropsi.f90
SRCS_0_HDF5   =	write_ogyropsi_hdf5.f90

# SRCS_1 chease_prog.f90 program wrapper needed if not in kepler
SRCS_1 = chease_prog.f90

# SRCS_2 get and put routines to read/write to ITM data structure, dummy routines if not needed
SRCS_2     = load_itm_dummy.f90 write_itm_dummy.f90
SRCS_2_ITM = load_itm_with_rout.f90 write_itm_with_rout.f90

# SRCS_interpos to make interpos library
# from: 
# svn export http://crppsvn.epfl.ch/repos/interpos/trunk/interpos_libs interpos_libs2
# cd interpos_libs2; cp -pr prec_rkind.f90 interpos_module.f90 cbsplgen.f90 cbsplgenp.f90 cbfitbnd.f90 cbfitper.f90 splipera.f90 splibnda.f90 intlinear.f90 intquadratic.f90 dpgbtrf_s.f90 ../interpos_libs
#
DIR_interpos = ./interpos_libs
MODS_interpos =	prec_rkind.f90 interpos_module.f90
SRCS_interpos =	cbsplgen.f90 cbsplgenp.f90 cbfitbnd.f90 cbfitper.f90 \
		splipera.f90 splibnda.f90 intlinear.f90 intquadratic.f90 
SRCS_interpos_lapack = dpgbtrf_s.f90
LIBINTERPOS=interpos

# For various targets:
# chease
SRCS		= $(SRCS_0) $(SRCS_0_NOHDF5) $(SRCS_1) $(SRCS_2)
# chease_HDF5
SRCS_HDF5	= $(SRCS_0) $(SRCS_0_HDF5) $(SRCS_1) $(SRCS_2)
# chease_ITM
SRCS_ITM	= $(SRCS_0) $(SRCS_0_NOHDF5) $(SRCS_1) $(SRCS_2_ITM)
# chease_ITM_HDF5
SRCS_ITM_HDF5	= $(SRCS_0) $(SRCS_0_HDF5) $(SRCS_1) $(SRCS_2_ITM)
# chease_KEPLER
SRCS_KEPLER	= $(SRCS_0) $(SRCS_0_NOHDF5) $(SRCS_2_ITM)

#INCS =	BNDIND.inc COMDAT.inc HERMIT.inc SOLOV.inc mdslib.inc
INCS =	BNDIND.inc COMDAT.inc HERMIT.inc SOLOV.inc

MODS     = euitm_schemas.f90 neobscoeffmod.f90 globals.f90 interpol.f90 prec_const.f90 \
	   sigmaneomod.f90 string_manipulation_tools.f90 euitm_xml_parser.f90
MODS_ITM = euitm_routines.f90
MODS_ITM = neobscoeffmod.f90 globals.f90 interpol.f90 prec_const.f90 \
	   sigmaneomod.f90 string_manipulation_tools.f90 euitm_xml_parser.f90

OBJS     = $(MODS:.f90=.o) $(SRCS:.f90=.o)
OBJS_HDF5     = $(MODS:.f90=.o) $(SRCS_HDF5:.f90=.o)
OBJS_ITM = $(MODS_ITM:.f90=.o) $(SRCS_ITM:.f90=.o)
OBJS_ITM_HDF5 = $(MODS_ITM:.f90=.o) $(SRCS_ITM_HDF5:.f90=.o)
OBJS_KEPLER = $(MODS:.f90=.o) $(MODS_ITM:.f90=.o) $(SRCS_KEPLER:.f90=.o)


.PRECIOUS:	$(SRCS) $(SRCS1)  $(SRCS1_ITM) $(INCS) $(MODS) $(MODS_ITM)

# Lutjens
# LIBS =	 -ldxml -lmat -lmx

# sun
LIBS =	 -xlic_lib=sunperf

# JET Linux
# LIBS =	 

# delphi
# LIBS =	 

CC = cc
CFLAGS = -O
FC = f77
FFLAGS = -O

F90 = f90

# delphi
#F90 = lf95

# F90FLAGS = -O0 -g
# LDFLAGS = -O0 -g

# Lutjens
#F90FLAGS = -fast -O -cpp -pipeline -arch ev6 -tune ev6
#LDFLAGS = -fast -O -cpp -pipeline -arch ev6 -tune ev6

# sun
F90FLAGS = -fast -O5
LDFLAGS = -fast -O5

# linux f95 delphi
#F90FLAGS = 
#LDFLAGS = -g

# Linux jac
# F90FLAGS = -fast -O4 -tp athlon
# LDFLAGS = -fast -O4 -tp athlon

# Nec SX7
F90 = sxf90
F90FLAGS =
LDFLAGS = 
LIBS = 

# Mac G5 IBM xlf90 compiler for Mac (Lutjens)
 F90 = xlf90
 F90FLAGS =  -qsuffix=f=f90 -O5 -qarch=g5 -qtune=g5 -qextname=mdsput2:mdsvalue:mdsconnect:mdsopen:mdsdisconnect
 LDFLAGS =
 LIBS = -L/usr/local/mdsplus/lib -lMdsLib

# Mac Intel ifort compiler for Mac (Lutjens)
 F90 = ifort
 F90FLAGS = -O3
 LDFLAGS =
 LIBS = -L/usr/local/mdsplus/lib -lMdsLib
 LIBS = 

# CRPPC74 AMD (Example using Futils and Hdf5)
# The Futils library is available on CRPP svn server
# 
#F90 = mpif90
#FUTILS = ${HOME}/PHD/Utils/hdf5/futils/src
#INCL_FUTILS = -I$(FUTILS)
#F90FLAGS = -O3
#HDF5=/usr/local/hdf5/lib
#INCL_HDF5 = -I${HDF5}
#LIBS=-L$(FUTILS) -lfutils \
#	-L$(HDF5) -lhdf5_fortran -lhdf5 -lz 

# g95 (on crpppc70)
F90 = g95
F90FLAGS = -r8  -ftrace=full -fno-second-underscore -fPIC 
LDFLAGS =
LIBS =	

# Note: if compiler does automatic deallocation euitm_xml_parser.f90 will not work. Specify here a lower optimization
F90FLAGS_parser = $(F90FLAGS)

# mpif90 (on enea142 gateway)
F90 = mpif90
F90FLAGS = -r8 -fPIC -Mnosecond_underscore -I/afs/efda-itm.eu/project/switm/ual/include/amd64_pgi 
LIBS = -L/afs/efda-itm.eu/project/switm/ual/lib/  -lUALFORTRANInterface_pgi  -L/afs/efda-itm.eu/project/switm/hdf5/amd64_pgi_1.8.1/lib/  -lz
F90FLAGS_parser = $(F90FLAGS)
INCL_HDF5 = 
INCL_FUTILS = -I /afs/efda-itm.eu/imp4/user/mcmillan/futils/1.2/src/

# g95 (on enea142 gateway)
F90 = g95
F90FLAGS = -r8  -ftrace=full -fno-second-underscore -fPIC 
F90FLAGS = -r8  -ftrace=full -fno-second-underscore -fPIC -I/afs/efda-itm.eu/project/switm/ual/include/amd64_g95
LDFLAGS =
LIBS = -L/afs/efda-itm.eu/project/switm/ual/lib/  -lUALFORTRANInterface_g95
F90FLAGS_parser = $(F90FLAGS)

# ifort (on crpppc70)
#F90 = ifort
#F90FLAGS = -c -O3 -r8 -automatic -xT -parallel
#LDFLAGS = -O3 -r8 -automatic -xT -parallel
#F90FLAGS_parser = -O1 -r8 -automatic -xT -parallel
#LIBS =	

# pleiades version
#F90 = mpif90
#FUTILS = ${HOME}/futils/src
#INCL_FUTILS = -I$(FUTILS)
#F90FLAGS = -O3
#F90FLAGS_parser = -O1
#HDF5=/usr/local/hdf5/lib
#INCL_HDF5 = -I${HDF5}
#LIBS=-L$(FUTILS) -lfutils \
#	-L$(HDF5) -lhdf5_fortran -lhdf5 -lz 

# Mac Intel ifort compiler for Mac (Lutjens)
F90 = ifort
F90FLAGS = -fast -xT -axT -i8
LDFLAGS =
LIBS = -L/usr/local/mdsplus/lib -lMdsLib
LIBS = -L/opt/intel/Compiler/11.0/064/Frameworks/mkl/lib/em64t -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_lapack -lmkl_core  -liomp5 -lguide -lm
F90FLAGS_parser = -fast -xT -axT -i8

# Nec SX8
F90 = sxf90
F90FLAGS = -Wf"-A dbl4" -C hopt -f2003
LDFLAGS = -Wf"-A dbl4" -C hopt -f2003
F90FLAGS = -ew -C hopt -f2003 
LDFLAGS = -ew -C hopt -f2003 
LIBS = 

# IBM (hal.epfl.ch) without HDF5 (at this stage do not compile file write_ogyropsi_hdf5.f90)
F90 = xlf90
F90FLAGS = -O3 -qtune=auto -qalign=4k -qdpc -q strict -qflag=E:E -qfree -qsuffix=f=f90
LDFLAGS = -O3 -qtune=auto
# for IBM need to compile without automatic deallocation of pointer, thus O<2 needed AND not xlf95
F90FLAGS_parser = -O1 -qtune=auto -qalign=4k -qdpc -q strict -qflag=E:E -qfree -qsuffix=f=f90
LIBS =	

F90 = mpif90
F90FLAGS = -r8 -fPIC -Mnosecond_underscore -I/afs/efda-itm.eu/project/switm/ual/include/amd64_pgi 
LIBS = -L/afs/efda-itm.eu/project/switm/ual/lib/  -lUALFORTRANInterface_pgi  -L/afs/efda-itm.eu/project/switm/hdf5/amd64_pgi_1.8.1/lib/  -lz
F90FLAGS_parser = $(F90FLAGS)
INCL_HDF5 = 
INCL_FUTILS = -I /afs/efda-itm.eu/imp4/user/mcmillan/futils/1.2/src/

# g95 (on enea142 gateway)
F90 = g95
F90FLAGS = -r8  -ftrace=full -fno-second-underscore -fPIC 
F90FLAGS = -r8  -ftrace=full -fno-second-underscore -fPIC -I/afs/efda-itm.eu/project/switm/ual/include/amd64_g95
LDFLAGS =
LIBS =  -Llib -L/afs/efda-itm.eu/project/switm/ual/lib/  -lUALFORTRANInterface_g95_ns
LIBS =  -Llib -L/afs/efda-itm.eu/project/switm/ual/lib/  -lUALFORTRANInterface_g95
F90FLAGS_parser = $(F90FLAGS)

# pgf90 (on enea142 gateway)
F90 = pgf90
F90FLAGS = -r8  -ftrace=full -fno-second-underscore -fPIC 
# use ueITM_schemas... from -I directory
F90FLAGS = -r8   -fPIC -I/afs/efda-itm.eu/project/switm/ual/4.06d/include/amd64_pgi
F90FLAGS = -r8   -fPIC -I/afs/efda-itm.eu/project/switm/ual/4.07a/include/amd64_pgi -I$(DIR_interpos)
LDFLAGS =
LIBS = -L/afs/efda-itm.eu/project/switm/ual/4.07a/lib  -lUALFORTRANInterface_pgi -L$(DIR_interpos) -l$(LIBINTERPOS)
F90FLAGS_parser = $(F90FLAGS)

# CSC desktop Warwick
F90 = mpif90
F90FLAGS = -r8 -O3 -I$(DIR_interpos)
LDFLAGS =
LIBS=-L$(DIR_interpos) -linterpos
F90FLAGS_parser = -O1

INCL_HDF5 = 
INCL_FUTILS = 

VERSION=V90_74


all: $(PROG)

$(PROG): $(DIR_interpos)/lib$(LIBINTERPOS).a $(OBJS)
# because of recursive call in euitm_xml_parser.f90 using pointers, 
# need compilation without automatic deallocation at this stage at least (June 2008)
	$(F90) $(F90FLAGS_parser) -c euitm_xml_parser.f90
	$(F90) $(F90FLAGS_parser) -c assign_code_parameters.f90
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

$(PROG_HDF5): $(DIR_interpos)/lib$(LIBINTERPOS).a $(OBJS_HDF5)
# because of recursive call in euitm_xml_parser.f90 using pointers, 
# need compilation without automatic deallocation at this stage at least (June 2008)
	$(F90) $(F90FLAGS_parser) -c euitm_xml_parser.f90
	$(F90) $(F90FLAGS_parser) -c assign_code_parameters.f90
	$(F90) $(LDFLAGS) -o $@ $(OBJS_HDF5) $(LIBS)

$(PROG_ITM): $(DIR_interpos)/lib$(LIBINTERPOS).a $(OBJS_ITM)
# so far chease_ITM can only be compiled on gateway, on which no need for special compilation for parser at this stage
	$(F90) $(LDFLAGS) -o $@ $(OBJS_ITM) $(LIBS)

$(PROG_ITM_HDF5): $(DIR_interpos)/lib$(LIBINTERPOS).a $(OBJS_ITM_HDF5)
# so far chease_ITM can only be compiled on gateway, on which no need for special compilation for parser at this stage
	$(F90) $(LDFLAGS) -o $@ $(OBJS_ITM_HDF5) $(LIBS)

$(PROG_KEPLER): $(DIR_interpos)/lib$(LIBINTERPOS).a $(OBJS_KEPLER)
# so far chease_ITM can only be compiled on gateway, on which no need for special compilation for parser at this stage
	ar -rv $(PROG_KEPLER).a $(OBJS_KEPLER)

$(PROG_interpos) :
	make -f Makefile_gateway $(DIR_interpos)/lib$(LIBINTERPOS).a

$(DIR_interpos)/lib$(LIBINTERPOS).a : 
	cd $(DIR_interpos); \
	$(F90) $(F90FLAGS)  -c $(MODS_interpos) $(SRCS_interpos) $(SRCS_interpos_lapack); \
	ar -rv lib$(LIBINTERPOS).a $(MODS_interpos:.f90=.o) $(SRCS_interpos:.f90=.o) $(SRCS_interpos_lapack:.f90=.o)
	pwd

segldr: 
	$(F90) $(LDFLAGS) -o chease $(OBJS) $(LIBS)

clean:
	rm -f  *.kmo *.mod *.M *.o

tar:
	gtar zcvf ../$(PROG).tar.gz $(SRCS) $(INCS) $(MODS) Makefile

source:	
	cat $(MODS) $(INCS) $(SRCS1) $(SRCS) > ../chease_$(VERSION).f90

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $< $(INCL_FUTILS) $(INCL_HDF5)

.f90.mod:
	@touch $*.mod
#	$(F90) $(F90FLAGS) -c $<

.mod.o:
	$(F90) $(F90FLAGS) -c $*.f90

acopy.o: globals.o
aldlt.o: prec_const.o
apcoef.o: globals.o
apcoef2.o: globals.o
assign_code_parameters.o: prec_const.o euitm_schemas.o euitm_xml_parser.o globals.o
#assign_code_parameters.o: prec_const.o euitm_xml_parser.o globals.o
atcoef.o: globals.o
auxval.o: globals.o
away.o: globals.o BNDIND.inc
ballit.o: globals.o
baloon.o: globals.o
basis1.o: globals.o HERMIT.inc
basis2.o: globals.o HERMIT.inc
basis3.o: globals.o HERMIT.inc
basis4.o: globals.o HERMIT.inc
blines.o: globals.o
bltest.o: globals.o interpol.o
bndspl.o: globals.o
bound.o: globals.o
bscoeff.o: globals.o neobscoeffmod.o
bsexpeq.o: globals.o
bsfunc.o: globals.o interpol.o
bstnzpro.o: globals.o
ccopy.o: globals.o
center.o: globals.o
check.o: globals.o BNDIND.inc
chease.o: globals.o
chipsi.o: globals.o
chipsimetrics.o: globals.o
cint.o: globals.o
conver.o: globals.o
copyap.o: globals.o interpol.o
copyapp.o: globals.o interpol.o
copyat.o: globals.o interpol.o
cotrol.o: globals.o
cubrt.o: globals.o
curent.o: globals.o
cvzero.o: globals.o
data.o: globals.o
direct.o: globals.o
drhodp.o: globals.o
dwy.o: prec_const.o
energy.o: globals.o
eqdim.o: globals.o
erdata.o: globals.o
errorch.o: globals.o
chease_prog.o: SOLOV.inc
euitm_schemas.o: SOLOV.inc
euitm_routines.o: euitm_schemas.o
euitm_xml_parser.o: prec_const.o string_manipulation_tools.o
evlate.o: globals.o
four1.o: prec_const.o
fourfft.o: globals.o
fourier.o: globals.o
g_0.o: globals.o
g_1.o: globals.o
g_2.o: globals.o
g_3.o: globals.o
gauss.o: globals.o
gchi.o: globals.o
gdataext.o: prec_const.o
genout.o: prec_const.o
gijlin.o: globals.o
gloadd.o: globals.o interpol.o
globals.o: prec_const.o
globals_init.o: globals.o
gloqua.o: globals.o interpol.o
guess.o: globals.o interpol.o
iarray.o: globals.o
identa.o: globals.o BNDIND.inc
identb.o: globals.o
indexx.o: prec_const.o
initia.o: globals.o
interpol.o: prec_const.o
iodisk.o: globals.o interpol.o COMDAT.inc
isamin.o: globals.o
ismax.o: globals.o
ismin.o: globals.o
isofind.o: globals.o
isofun.o: globals.o
isrchfge.o: globals.o
issum.o: globals.o
itipr.o: globals.o
ivar.o: globals.o
jnovaw.o: globals.o
labrun.o: globals.o
limita.o: globals.o
limitb.o: globals.o
load_itm_dummy.o: globals.o
load_itm_with_rout.o: globals.o
ltxw.o: prec_const.o
lyv.o: prec_const.o
magaxe.o: globals.o
mappin.o: globals.o interpol.o
matrix.o: globals.o
mesage.o: globals.o
mesh.o: globals.o
msplcy.o: globals.o
mspline.o: globals.o
nerat.o: globals.o interpol.o
nonlin.o: globals.o
norept.o: globals.o interpol.o
ntridg.o: globals.o
oarray.o: globals.o
oldeq.o: globals.o
oldnew.o: globals.o
outgload.o: globals.o interpol.o
outmetric.o: globals.o
outmksa.o: globals.o
outnvw.o: globals.o interpol.o
outpen.o: globals.o
output.o: globals.o interpol.o COMDAT.inc
outxt.o: globals.o
packme.o: globals.o interpol.o
packmep.o: globals.o
page.o: globals.o
polyfun.o: globals.o
polynm.o: globals.o interpol.o
ppbstr.o: globals.o interpol.o
pprime.o: globals.o interpol.o
pprm.o: globals.o
ppspln.o: globals.o prec_const.o
ppspln2.o: prec_const.o
premap.o: globals.o interpol.o
preset.o: globals.o
prfunc.o: globals.o
priqqu.o: globals.o
prnorm.o: globals.o
profile.o: globals.o interpol.o
psibox.o: globals.o interpol.o
psicel.o: globals.o
psvol.o: globals.o interpol.o
qplacs.o: globals.o interpol.o
rarray.o: globals.o
realft.o: prec_const.o
reseti.o: globals.o
resetr.o: globals.o
resppr.o: globals.o
rmrad.o: globals.o
rscale.o: globals.o
rvar.o: globals.o
rvar2.o: globals.o
rzbound.o: globals.o
runtim.o: prec_const.o
saxpy.o: prec_const.o
scopy.o: prec_const.o
scopyr.o: prec_const.o
sdot.o: prec_const.o
setupa.o: globals.o
setupb.o: globals.o
shave.o: globals.o
smooth.o: globals.o
solovev.o: globals.o SOLOV.inc
solvit.o: globals.o
sort3.o: prec_const.o
splcy.o: globals.o
splcyp.o: globals.o
splifft.o: prec_const.o
spline.o: prec_const.o interpol.o
sscal.o: prec_const.o
ssum.o: globals.o
stchps.o: globals.o
stepon.o: globals.o interpol.o
string_manipulation_tools.o: prec_const.o
subsz.o: globals.o
surface.o: globals.o interpol.o sigmaneomod.o
surfadd.o: globals.o interpol.o
surfrz.o: globals.o
tcase.o: globals.o
test.o: globals.o SOLOV.inc
tetare.o: globals.o
tpsi.o: globals.o
tricyc.o: prec_const.o
tricycm.o: prec_const.o
tridagm.o: prec_const.o
tshift.o: globals.o
vacufft.o: globals.o
vacuum.o: globals.o
vlion.o: globals.o
vzero.o: prec_const.o
write_itm_dummy.o: globals.o
write_itm_with_rout.o: globals.o
wrtext.o: globals.o
wrtplot.o: globals.o interpol.o COMDAT.inc
wrtbin.o: globals.o interpol.o COMDAT.inc
wrtasc.o: globals.o interpol.o COMDAT.inc
wrtmat.o: globals.o interpol.o COMDAT.inc
wrtbin.o: globals.o interpol.o COMDAT.inc
outgyro.o: globals.o interpol.o
outelit.o: globals.o
outastro.o: globals.o interpol.o
outgyropsi.o: globals.o
outgyropsi_hdf5.o:globals.o
hamada.o: globals.o interpol.o
neoart.o: globals.o
