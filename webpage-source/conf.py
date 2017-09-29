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
source_suffix = ['.rst', '.md']
source_parsers = {
   '.md': 'recommonmark.parser.CommonMarkParser',
}
master_doc = 'index'

# General information about the project.
project = 'acados'
copyright = '2017, Syscop'
author = 'Dimitris Kouzoupis, Robin Verschueren, Gianluca Frison'

version = '0.5'
release = '0.5a1'
language = None

exclude_patterns = ['_build', 'README.md', 'Thumbs.db', '.DS_Store']
pygments_style = 'sphinx'
todo_include_todos = True

# -- Options for HTML output ----------------------------------------------
html_theme = 'alabaster'
html_theme_options = {
    'github_user': 'acados',
    'github_repo': 'acados',
    'github_banner': True
}
html_sidebars = {
    '**': [
        'about.html',
        'navigation.html',
        'relations.html'  # needs 'show_related': True theme option to display
    ]
}
html_static_path = ['webpage-static']
html_copy_source = False
html_sourcelink_suffix = ''
