SOURCES = $(wildcard *.Rnw)
FINALS = $(patsubst %.Rnw, %.pdf, $(SOURCES))
$(info Running on $(HOSTNAME))

all: $(FINALS)

%.pdf: %.asy
	asy -f pdf $<

intro-to-maximum-likelihood.pdf: probability-density.pdf

%.pdf: %.Rnw
	echo $HOSTNAME && Sweave.sh -ld $<

.PHONY: all
