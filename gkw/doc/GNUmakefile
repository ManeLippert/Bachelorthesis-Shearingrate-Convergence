.PHONY=all doc doc-diff devel derivation intro ballooning gs2gkw gkwandchease runchease source clean
# where to look for .tex files
#VPATH = ./:./etc/:./manual/:./chease/
vpath %.tex ./etc/ ./manual/ ./chease/
vpath %.sty ./etc/ ./manual/ ./chease/
vpath %.dvi ./etc/ ./manual/ ./chease/
vpath %.ps ./etc/ ./manual/ ./chease/

#BUILD_FOLDER=./build

LATEX = latex -src
DVIPS = dvips -Ppdf -G0
PS2PDF = ps2pdf14 -sPAPERSIZE=a4
DVIPDFM = dvipdfm -p a4
PDFLATEX = pdflatex
RERUN_LATEX_LINE='may have changed. Rerun to get\|Rerun to get outlines right'
FIGURES = $(shell find fig \( -name '*eps' -o -name '*tex' \) -type f)

all: doc devel derivation intro ballooning gs2gkw gkwandchease runchease

# $(BUILD_FOLDER):
# 	mkdir -p $@
# 	mkdir -p $@/manual

doc: manual/doc.pdf $(FIGURES) $(shell find benchmarks -name '*eps' -type f)

doc-diff:
	@echo "Note: The doc-diff target does not see any uncommitted text. You probably want to (temporarily) commit everything locally."
	@echo "The doc-diff target cannot compare to commits earlier than sha 193501014'"
	@read -p "Compare manual to which git commit? Enter sha,tag,branch or similar: " REPLY ;\
	cd manual ; \
	git latexdiff  --output doc-diff.pdf --ignore-latex-errors --verbose --flatten --latex --main doc.tex $$REPLY --
	@echo
	@echo "Please see the file manual/doc-diff.pdf !"


devel: etc/devel.pdf

derivation: etc/derivation.pdf  $(FIGURES)

intro: etc/GKW_for_beginners.pdf $(FIGURES)

ballooning: etc/Ballooning.pdf

gs2gkw: etc/gs2gkw.pdf

gkwandchease: chease/gkwandchease.pdf

runchease: chease/runchease.pdf

source:
	@echo "Start extracting comments from source files and making the html/pdf output."
	@echo " Output of doxygen/latex is written to corresponding log/err file."
	@echo "---------------------------------------------------------------------------"
	doxygen 1> doxygen.log 2> doxygen.err
	make -C latex/ 1> latex.log 2> latex.err
	cp latex/refman.pdf .
	make -C latex/ clean
	@echo "---------------------------------------------------------------------------"
	@echo "...done."

clean:
	rm -f {.,etc,chease}/*.pdf
	rm -f manual/{doc,doc-diff}.pdf
	rm -f {.,manual,etc,chease}/*.{ps,log,aux,nlo,toc,dvi,out}
	rm -rf latex html


################################################################################
### rules                                                                    ###
################################################################################

%.dvi : %.tex
	@echo make $@ with $< in $(<D) ;\
	cd $(<D); \
	$(LATEX) $(<F) ; \
	while [ -n "$$(grep $(RERUN_LATEX_LINE) *.log)" ]; do \
	     $(LATEX) $(<F); \
	done

	# if [ -e $(BUILD_FOLDER)/$*.nlo -a -e  $(BUILD_FOLDER)/nomencl.ist ]; then \
	#   makeindex $(BUILD_FOLDER)/$*.nlo -s $(BUILD_FOLDER)/nomencl.ist -o $(BUILD_FOLDER)/$*.nls ; \
	# fi
	# if [ -e $(BUILD_FOLDER)/$*.idx ]; then \
	#   makeindex $(BUILD_FOLDER)/$*.idx ; \
	# fi 

%.ps : %.dvi
	@echo make $@ with $< in $(<D) ;\
	cd $(<D) ; pwd ;\
	$(DVIPS) $(<F) -o $(@F)

%.pdf : %.ps
	@echo make $@ with $< in $(<D);\
	cd $(<D) ; pwd ;\
	$(PS2PDF) $(<F) $(@F)

#%.pdf : %.tex
#	$(PDFLATEX) $<


manual/doc.pdf : preface.tex theory.tex practise.tex collisions.tex rotation.tex \
          tearing.tex implementation.tex buildandrun.tex diagnostics.tex \
          benchmarks.tex bibliography.tex Appendixgeometry.tex eiv.tex \
          Memo-LFS-FSA-density.tex global.tex gkwdoccommands.sty
