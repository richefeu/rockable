# Build the technical note and copy the PDF where Sphinx offers it for download.
MAIN = skin-meshes

.PHONY: all clean

all: ../$(MAIN).pdf

../$(MAIN).pdf: $(MAIN).tex $(MAIN).bib $(wildcard figures/*)
	latexmk -pdf -halt-on-error -interaction=nonstopmode $(MAIN).tex
	cp $(MAIN).pdf ../$(MAIN).pdf

clean:
	latexmk -C $(MAIN).tex
	rm -f $(MAIN).bbl
