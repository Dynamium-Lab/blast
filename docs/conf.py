# Configuration file for the Sphinx documentation builder.
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

# -- Project information -----------------------------------------------------

project = "Blast"
copyright = "2024-2026, Andre Gallant"
author = "Andre Gallant"
release = "1.0"

# -- General configuration ---------------------------------------------------

extensions = [
    "myst_parser",          # author pages in Markdown
    "breathe",              # C++ API via Doxygen XML
    "sphinx.ext.autodoc",   # Python API from docstrings
    "sphinx.ext.napoleon",  # NumPy/Google-style docstring parsing
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "_doxygen", "Thumbs.db", ".DS_Store"]

# MyST: allow ::: fenced directives and useful extensions.
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "fieldlist",
]
myst_heading_anchors = 3

# -- Breathe (C++) -----------------------------------------------------------

breathe_projects = {"blast": "_doxygen/xml"}
breathe_default_project = "blast"
breathe_default_members = ("members", "undoc-members")

# -- Autodoc (Python) --------------------------------------------------------
# The compiled `blast` extension must be importable at build time
# (RTD installs it via pip; locally run `pip install .` first).

autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
}
autodoc_member_order = "groupwise"
napoleon_numpy_docstring = True
napoleon_google_docstring = True

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
}

# -- HTML output -------------------------------------------------------------

html_theme = "furo"
html_title = "Blast"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_show_sphinx = False

html_theme_options = {
    "source_repository": "https://github.com/Dynamium-Lab/blast/",
    "source_branch": "main",
    "source_directory": "docs/",
    "light_css_variables": {
        "color-brand-primary": "#3b5bdb",
        "color-brand-content": "#3b5bdb",
    },
    "dark_css_variables": {
        "color-brand-primary": "#aebbfe",
        "color-brand-content": "#aebbfe",
    },
}
