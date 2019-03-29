# Acados Sphinx Makefile
# --
# Based on the sphinx-quickstart Makefile template.
# Author: Rezart Qelibari <qelibarr@informatik.uni-freiburg.de>

SPHINXOPTS    =
SPHINXBUILD   = python -msphinx
SPHINXPROJ    = acados
SOURCEDIR     = .
BUILDDIR      = _build
GENERATED     = .doctrees _static _sources .buildinfo .nojekyll genindex.html objects.inv search.html searchindex.js

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile clean html

clean:
	@for file in *.rst; do rm "$(BUILDDIR)/$${file%.*}.html" 2> /dev/null || true; done
	@for file in $(GENERATED); do CPATH="$(BUILDDIR)/$${file}"; [ -f "$${CPATH}" -o -d "$${CPATH}" ] && rm -r "$${CPATH}" || true; done

html:
	@$(SPHINXBUILD) -b $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
