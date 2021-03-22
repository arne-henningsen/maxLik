SOURCES = $(wildcard *.Rnw)
FINALS = $(patsubst %.Rnw, %.pdf, $(SOURCES))

all: $(FINALS)

%.pdf: %.asy
	asy -f pdf $<

intro-to-maximum-likelihood.pdf: probability-density.pdf

%.pdf: %.Rnw
	Sweave.sh -ld $<

.PHONY: all
