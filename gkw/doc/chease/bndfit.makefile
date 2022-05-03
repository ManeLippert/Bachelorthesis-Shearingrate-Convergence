SHELL=/bin/sh
######################################################################
#  makefile           created on Tue Mar 14 17:21:34 MET 1995
#
#   This file is an utility to compile the program bndfit. It has
# been generated using the program script mkgen V1.5.2.
#
# It is envoked by the command line
#	make bndfit
#  If the source preprocessor was selected, the GENS files are first 
# transformed into fortran source by changing lines which start with
#       C+ FORTRAN COMMAND
# into their active equivalent (without C+), following the 
# rules of the awk preprocessor .
#
#  Compilation of the fortran modules referred to as FORS into the object
# modules OBJS and then loading to produce the executable bndfit is 
# taken care of by cf77 -c. The resulting file is copied to the 
# directory containing also the make directory MK.
#
# The following COMMANDs may be constructed:
#   make bndfit	: preprocess & compile modules which are not up to
#                       date and create executable
#   make source		: to create a source file of the code with include files
#   make unsplit	: recreates the source program by sticking all 
#                       the modules together and delete subdir MK
#   make listing        : makes a listing in bndfit.l including line
#                         number for each module and for total file
#                        for partial listing, change  in makefile
#                         (assume pagenumberlist.exe known command,
#                          should in same place as mkgen)
#
# Cray2 allows to select the following compiler options:
#   make bndfit.flow	: compile & load with flowtrace options
#   make bndfit.prof	: compile & load with profview options
#   make bndfit.cdbx	: compile & load with interactive debugging opt
#
######################################################################
#
# --- Tunable parameters ---------------------------------------------
#
PREPROC =
MACHINE =
#
BUILDER =f77
COMPILER=f77 -c
#COMPOPT = -Wf" -ez -ooff"
COMPOPT = -Wf" -ez"
#
# sun and sgi
COMPOPT = -r8 -i4
#
LOADOPT = -r8 -i4
LOADOPT = 

# sun
LOADLIB = -xlic_lib=sunperf
# sgi
#LOADLIB = -lcomplib.sgimath
# Linux-jac Athlon
LOADLIB = 
#
COMMAND  =bndfit
COMMAND_eff  =bndfit
# Version number, example: VERSION = .A940101
VERSION = .B980526
#
# --- Define rules ---------------------------------------------------
#
#.SUFFIXES :	
#.SUFFIXES :	.f .f
#
#  --- Preprocessor
#.f.f :	
#	- rm $*.f
#	awk -f $(PREPROC) machine=$(MACHINE) $<   >$*.f
#
#  --- Compiler
.f.o :	
	$(COMPILER) $(COMPOPT)  $<
#
#
#  --- Macro definitions ---------------------------------------------
#

INCS=		 \
	QUAQQQ.inc

FORS=	main.f  bndinp.f  equimsh.f  interp.f  spline.f  splint.f   \
	polfit.f  gauss.f  fftfit.f  sincosft.f  polout.f  sort3.f  indexx.f  ismax.f  ismin.f   \
	mulplt.f  mulpld.f  mulpdo.f  graphe.f  graphn.f  graphnm.f  gradot.f   \
	grapdn.f  grallg.f  grllgn.f  gralgl.f  grlglg.f  glashn.f  gcushn.f   \
	gcusdn.f  gxigrn.f  gxishn.f  pgtxnb.f  gtx.f  gchar.f  gvect.f   \
	gdot.f  gclrwk.f  plots.f  plotf.f   bndanal.f cbsplgen_per.f cubfit.f nonsym.f

FORLIS=	main.f  bndinp.f  equimsh.f  interp.f  spline.f  splint.f   \
	polfit.f  gauss.f  fftfit.f  sincosft.f  polout.f  sort3.f  indexx.f  ismax.f  ismin.f   \
	mulplt.f  mulpld.f  mulpdo.f  graphe.f  graphn.f  graphnm.f  gradot.f   \
	grapdn.f  grallg.f  grllgn.f  gralgl.f  grlglg.f  glashn.f  gcushn.f   \
	gcusdn.f  gxigrn.f  gxishn.f  pgtxnb.f  gtx.f  gchar.f  gvect.f   \
	gdot.f  gclrwk.f  plots.f  plotf.f  bndanal.f cbsplgen_per.f cubfit.f nonsym.f

OBJS=	main.o  bndinp.o  equimsh.o  interp.o  spline.o  splint.o   \
	polfit.o  gauss.o fftfit.o  sincosft.o polout.o  sort3.o  indexx.o  ismax.o  ismin.o   \
	mulplt.o  mulpld.o  mulpdo.o  graphe.o  graphn.o  graphnm.o  gradot.o   \
	grapdn.o  grallg.o  grllgn.o  gralgl.o  grlglg.o  glashn.o  gcushn.o   \
	gcusdn.o  gxigrn.o  gxishn.o  pgtxnb.o  gtx.o  gchar.o  gvect.o   \
	gdot.o  gclrwk.o  plots.o  plotf.o  bndanal.o cbsplgen_per.o cubfit.o nonsym.o

#
#  -- Files not to be deleted when make is interrupted
#
.PRECIOUS:	$(FORS) $(INCS)
#
#
# --- Commands -------------------------------------------------------
#
# Lines from here on down should not need to be changed.  They are the
# actual rules which make uses to build bndfit in different ways.
#
#
# --- Executable (first command so is executed if only "make" is typed)
#
$(COMMAND)	:  $(OBJS)
	$(BUILDER) $(LOADOPT) -o $(COMMAND) $(OBJS) $(LOADLIB)
#	touch $(COMMAND) $(COMMAND).flow $(COMMAND).prof $(COMMAND).dbg
#	mv $(COMMAND) ../$(COMMAND_eff)
#	touch $(COMMAND)
#
# --- touch .f file if a referenced include file is younger
#
$(FORS)	: $(INCS)
	@for x in $?; do \
	  if [ `grep $$x $@ | wc -l` -ne 0 ]; then \
	    touch $@; \
	  fi \
	done
# --- Common blocks
#
#$(FORS)	: $(INCS)
#	touch $@
#
# --- Flowtrace (not available on the SUN)
#
$(COMMAND).flow :  $(FORS)
	$(BUILDER) -c -F -Wf"-em" $(COMPOPT) `echo " $?" `
	$(BUILDER) $(LOADOPT) -o $(COMMAND).flow   $(OBJS) $(LOADLIB)
	touch $(COMMAND) $(COMMAND).flow $(COMMAND).prof $(COMMAND).dbg
	mv $(COMMAND).flow ../$(COMMAND).flow
	touch $(COMMAND).flow
#
# --- Profview
#
$(COMMAND).prof :  $(FORS)
	$(BUILDER) -c -Wf"-emz" $(COMPOPT) `echo " $?" `
	$(BUILDER) $(LOADOPT) -o $(COMMAND).prof -lprof $(LOADLIB) $(OBJS)
	touch $(COMMAND) $(COMMAND).flow $(COMMAND).prof $(COMMAND).dbg
	mv $(COMMAND).prof ../$(COMMAND).prof
	touch $(COMMAND).prof
#
# --- Debbugger (symbol table + check array boundaries)
#
$(COMMAND).dbg	:  $(FORS)
	$(BUILDER) -c -Wf"-eimz -Rbc" $(COMPOPT) `echo " $?" `
	$(BUILDER) $(LOADOPT) -o $(COMMAND).dbg   $(OBJS) $(LOADLIB)
	touch $(COMMAND) $(COMMAND).prof $(COMMAND).dbg
	mv $(COMMAND).dbg ../$(COMMAND).dbg
	touch $(COMMAND).dbg
#
source	:
	cat $(INCS) $(FORS) >../$(COMMAND)$(VERSION).src
#
unsplit	:
	cat $(FORS) >../$(COMMAND)$(VERSION).f
	- cp $(INCS) ..
	################################################################
	#===> You may now remove safely the makefile directory by typing
	#        cd .. ; rm -r MK                                   
	################################################################

#
tar	:
	@mkdir MK$(VERSION)
	@mv *.inc *.f  makefile MK$(VERSION)
	tar cf ../$(COMMAND)$(VERSION)tar MK$(VERSION)
	@mv MK$(VERSION)/* .
	@rm -r MK$(VERSION)

listing	:
#	/bin/cc -opagenumberlist.exe pagenumberlist.c
	@echo " " > sou.update
	@rm sou.update
	@for x in $(INCS); do \
	  echo "@COMDECK "$$x >> sou.update; \
	  cat $$x >> sou.update; \
	done
	@for x in $(FORLIS); do \
	  echo "@DECK "$$x  >> sou.update; \
	  cat $$x >> sou.update; \
	done
#     assume pagenumberlist.exe known command (should in same place as mkgen)
	@echo "********  create listing ***************"
	pagenumberlist.exe < sou.update >../$(COMMAND)$(VERSION).src.l
	@rm sou.update
