SOURCES = $(wildcard *.Rnw)
PDFS = $(patsubst %.Rnw, %.pdf, $(SOURCES))
FINALS = $(patsubst %.pdf, ../inst/doc/%.pdf, $(PDFS))
SWEAVE = /usr/lib/R/bin/Sweave

.PHONY: all
all: $(PDFS) $(FINALS)

%.pdf: %.asy
	asy -f pdf $<

intro-to-maximum-likelihood.pdf: probability-density.pdf

%.tex: %.Rnw
	$(SWEAVE) --clean $<

%.pdf: %.tex
	latexmk -xelatex -f -interaction=nonstopmode $<

../inst/doc/%.pdf: %.pdf
	cp $< ../inst/doc/

.PHONY: clean
clean:
	rm -vf *.tex *.aux *.log *.R *.blg *.fdb_latexmk\
	 *.fls *.xdv *.bbl *.out
