# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build
PYTHON        ?= python

# The version that versioned-html builds; the publish workflow sets it to
# develop or to the release version
DOCS_VERSION  ?= develop
export DOCS_VERSION

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile versioned-html versioned-html-strict \
	versioned-shared clean-versioned-html

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# raise warnings to errors
html-strict:
	@$(SPHINXBUILD) -b html -nW --keep-going "$(SOURCEDIR)" "$(BUILDDIR)/html" $(SPHINXOPTS) $(O)

# Build into _build/html/$(DOCS_VERSION) with the version switcher, laid out
# as on GitHub Pages
versioned-html:
	@$(SPHINXBUILD) -b html "$(SOURCEDIR)" "$(BUILDDIR)/html/$(DOCS_VERSION)" $(SPHINXOPTS) $(O)
	@$(MAKE) --no-print-directory versioned-shared

# The same, raising warnings to errors
versioned-html-strict:
	@$(SPHINXBUILD) -b html -nW --keep-going "$(SOURCEDIR)" "$(BUILDDIR)/html/$(DOCS_VERSION)" $(SPHINXOPTS) $(O)
	@$(MAKE) --no-print-directory versioned-shared

# The switcher script and the list of versions built so far
versioned-shared:
	mkdir -p "$(BUILDDIR)/html/shared"
	cp shared/version-switcher.js "$(BUILDDIR)/html/shared/"
	$(PYTHON) generate_versions_json.py --base-dir "$(BUILDDIR)/html"
	@echo "To preview, run 'python -m http.server -d $(BUILDDIR)/html' and open"
	@echo "http://localhost:8000/$(DOCS_VERSION)/index.html"

clean-versioned-html:
	rm -rf $(BUILDDIR)/html/*

clean:
	rm -rf $(BUILDDIR)
