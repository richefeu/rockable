# Build the Rockable cheatsheets (Typst).
#
#   make            all the PDFs
#   make watch SHEET=physics_rockable_cheatsheet
#   make clean

TYPST ?= typst

SOURCES := $(filter-out common.typ,$(wildcard *.typ))
PDFS    := $(SOURCES:.typ=.pdf)

SHEET ?= configuration_rockable_cheatsheet

.PHONY: all watch clean

all: $(PDFS)

%.pdf: %.typ common.typ
	$(TYPST) compile $<

watch: $(SHEET).typ common.typ
	$(TYPST) watch $(SHEET).typ

clean:
	rm -f $(PDFS)
