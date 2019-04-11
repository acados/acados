# Acados Sphinx Makefile
# --
# Based on the sphinx-quickstart Makefile template.

SPHINXOPTS    =
SPHINXBUILD   = python -msphinx
DOXYBUILD     = doxygen
SPHINXPROJ    = acados
SOURCEDIR     = .
BUILDDIR      = _build
GENERATED     = .doctrees _static _sources .buildinfo .nojekyll genindex.html objects.inv search.html searchindex.js
DOXYBUILDDIR  = _build_doxygen_c_api _build_doxygen


default: doxygen html

help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

html:
	@$(SPHINXBUILD) -b $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

doxygen:
	@$(DOXYBUILD) c_api/Doxyfile
	@$(DOXYBUILD) doxygen/Doxyfile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

clean:
	@echo Clean Sphinx
	@for file in *.rst; do rm "$(BUILDDIR)/$${file%.*}.html" 2> /dev/null || true; done
	@for file in $(GENERATED); do CPATH="$(BUILDDIR)/$${file}"; [ -f "$${CPATH}" -o -d "$${CPATH}" ] && rm -r "$${CPATH}" || true; done
	@echo Clean doxygen
	@rm -r $(DOXYBUILDDIR)

.PHONY: help Makefile clean html default doxygen
