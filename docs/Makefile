#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

# Acados Sphinx Makefile
# --
# Based on the sphinx-quickstart Makefile template.

SPHINXOPTS    =
SPHINXBUILD   = python -msphinx
SPHINXBUILD   = sphinx-build
DOXYBUILD     = doxygen
SPHINXPROJ    = acados
SOURCEDIR     = .
BUILDDIR      = _build
GENERATED     = .doctrees _static _sources .buildinfo .nojekyll genindex.html objects.inv search.html searchindex.js
DOXYBUILDDIR  = _build_doxygen # _build_doxygen_c_interface


default: doxygen html

help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

html:
	@$(SPHINXBUILD) -b $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

doxygen:
	@$(DOXYBUILD) doxygen/Doxyfile
	@$(DOXYBUILD) c_interface/Doxyfile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

upload:
	rsync -aP -e ssh $(BUILDDIR)/ acados_website@syscop.de:~/website_public/ --delete

clean:
	@echo Clean Sphinx
	@for file in *.rst; do rm "$(BUILDDIR)/$${file%.*}.html" 2> /dev/null || true; done
	@for file in $(GENERATED); do CPATH="$(BUILDDIR)/$${file}"; [ -f "$${CPATH}" -o -d "$${CPATH}" ] && rm -r "$${CPATH}" || true; done
	@echo Clean doxygen
	@rm -r $(DOXYBUILDDIR)
	@rm -r $(BUILDDIR)


.PHONY: help Makefile clean html default doxygen
