# Makefile for Sphinx documentation
#

SHELL := /bin/bash

# You can set these variables from the command line.
SPHINXOPTS    = -v
SPHINXBUILD   ?= LANG=C sphinx-build
PAPER         =

BUILD_DIR     ?= build
FILES         =

# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -d $(BUILD_DIR)/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) source

.PHONY: help clean html

#------------------------------------------------------------------------------
all: html
help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html      to make standalone HTML files"

clean:
	-rm -rf $(BUILD_DIR)/* source/_generated source/api/_generated

html:
	mkdir -p $(BUILD_DIR)/html $(BUILD_DIR)/doctrees
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILD_DIR)/html $(FILES)
	@echo
	@echo "Build finished. The HTML pages are in build/html."

html_check:
	SIRE_NO_LAZY_IMPORT=1 $(SPHINXBUILD) -b linkcheck $(ALLSPHINXOPTS) build/html $(FILES)

