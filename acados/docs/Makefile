#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# The 2-Clause BSD License
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.;
#

# Acados Sphinx Makefile
# --
# Based on the sphinx-quickstart Makefile template.

SPHINXOPTS    =
# SPHINXBUILD   = python -msphinx
SPHINXBUILD   = sphinx-build
DOXYBUILD     = doxygen
SPHINXPROJ    = acados
SOURCEDIR     = .
BUILDDIR      = _build
GENERATED     = .doctrees _static _sources .buildinfo .nojekyll genindex.html objects.inv search.html searchindex.js
DOXYBUILDDIR  = _build_doxygen_c_interface _build_doxygen


default: doxygen html

help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

html:
	@$(SPHINXBUILD) -b $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

doxygen:
	@$(DOXYBUILD) c_interface/Doxyfile
	@$(DOXYBUILD) doxygen/Doxyfile

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
