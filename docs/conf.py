# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import subprocess

import recommonmark
from recommonmark.transform import AutoStructify

source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'markdown',
    '.md': 'markdown',
}

templates_path = ['_templates']
master_doc = 'index'

breathe_projects = { "acados": "_build_doxygen_c_api/xml/" }
breathe_default_project = "acados"
subprocess.call('doxygen c_api/Doxyfile', shell=True)
subprocess.call('doxygen doxygen/Doxyfile', shell=True)

# -- Project information -----------------------------------------------------

project = 'acados'
copyright = '2019, syscop'
author = 'syscop'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [ 'breathe', 'recommonmark' ]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

exclude_patterns = ['_build', 'README.md', 'Thumbs.db', '.DS_Store']
pygments_style = 'sphinx'
todo_include_todos = True

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'
# html_theme = 'alabaster'

html_theme_options = {
    'canonical_url': '',
    'analytics_id': 'UA-XXXXXXX-1',  #  Provided by Google in your dashboard
    'logo_only': False,
    'display_version': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    # 'vcs_pageview_mode': '',
    # Toc options
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False
}
# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

def setup(app):
    app.add_config_value('recommonmark_config', {
            'enable_auto_toc_tree': True,
            'auto_toc_tree_section': 'Contents',
            'enable_eval_rst': True,
            'enable_math':True
            }, True)
    app.add_transform(AutoStructify)
