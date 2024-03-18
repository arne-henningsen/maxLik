SOURCES = $(wildcard *.Rnw)
FINALS = $(patsubst %.Rnw, %.pdf, $(SOURCES))
SWEAVE = /usr/lib/R/bin/Sweave
$(info Running on $(HOSTNAME))

all: $(FINALS)

%.pdf: %.asy
	asy -f pdf $<

intro-to-maximum-likelihood.pdf: probability-density.pdf

%.pdf: %.Rnw
	$(SWEAVE) --pdf --clean $<

.PHONY: all
