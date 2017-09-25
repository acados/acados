#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Acados Sphinx Config
# --
# Author: Rezart Qelibari <qelibarr@informatik.uni-freiburg.de>

# -- General configuration ------------------------------------------------
extensions = ['sphinx.ext.todo',
    'sphinx.ext.mathjax',
    'sphinx.ext.githubpages']
templates_path = ['webpage-templates']
source_suffix = '.rst'
master_doc = 'index'

# General information about the project.
project = 'acados'
copyright = '2017, D. Kouzoupis et al'
author = 'Dimitris Kouzoupis, Robin Verschueren, Gianluca Frison'

version = '0.5'
release = '0.5a1'
language = None

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
pygments_style = 'sphinx'
todo_include_todos = True

# -- Options for HTML output ----------------------------------------------
html_theme = 'alabaster'
html_sidebars = {
    '**': [
        'about.html',
        'navigation.html',
        'relations.html',  # needs 'show_related': True theme option to display
        'searchbox.html',
        'donate.html'
    ]
}
html_static_path = ['webpage-static']
