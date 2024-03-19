SOURCES = $(wildcard *.Rnw)
PDFS = $(patsubst %.Rnw, %.pdf, $(SOURCES))
FINALS = $(patsubst %.pdf, ../inst/doc/%.pdf, $(PDFS))
SWEAVE = /usr/lib/R/bin/Sweave

.PHONY: all
all: $(FINALS)

%.pdf: %.asy
	asy -f pdf $<

intro-to-maximum-likelihood.pdf: probability-density.pdf

%.pdf: %.Rnw
	$(SWEAVE) --pdf --clean $<

../inst/doc/%.pdf: %.pdf
	mv $< ../inst/doc/

.PHONY: clean
clean:
	rm -vf *.tex *.aux *.log *.R *.blg *.fdb_latexmk\
	 *.fls *.xdv *.bbl *.out *.pdf
