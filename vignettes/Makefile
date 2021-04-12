SOURCES = $(wildcard *.Rnw)
FINALS = $(patsubst %.Rnw, %.pdf, $(SOURCES))
$(info Running on $(HOSTNAME))

all: $(FINALS)

%.pdf: %.asy
	asy -f pdf $<

intro-to-maximum-likelihood.pdf: probability-density.pdf

# on r-forge Sweave.sh does not work...
ifneq ($(HOSTNAME),)
%.pdf: %.Rnw
	Sweave.sh -ld $<
endif

.PHONY: all
