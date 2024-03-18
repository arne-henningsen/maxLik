SOURCES = $(wildcard *.Rnw)
PDFS = $(patsubst %.Rnw, %.pdf, $(SOURCES))
FINALS = $(patsubst %.pdf, ../inst/doc/%.pdf, $(PDFS))
SWEAVE = /usr/lib/R/bin/Sweave
$(info Running on $(HOSTNAME))

all: $(FINALS)

%.pdf: %.asy
	asy -f pdf $<

intro-to-maximum-likelihood.pdf: probability-density.pdf

%.pdf: %.Rnw
	$(SWEAVE) --pdf --clean $<

../inst/doc/%.pdf: %.pdf
	mv $< ../inst/doc/

.PHONY: all
