# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Blast'
copyright = '2024, Andre Gallant'
author = 'Andre Gallant'
release = '1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['myst_parser']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
html_show_sphinx = False
html_show_sourcelink = False
html_copy_source = False

pygments_style = 'nord'


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_theme_options = {
    'fixed_sidebar': 'True',
    'page_width': '1100px',
    'description' : 'A fast C++ trajectory optimization library for robot manipulators',
    'base_bg': '#1E2129',
    'gray_1' : '#BEC1C6',
    'base_text' : '#BEC1C6',
    'body_text' : '#BEC1C6',
    'footer_text' : '#BEC1C6',
    'sidebar_text' : '#BEC1C6',
    'sidebar_list' : '#BEC1C6',
    'sidebar_hr' : '#000',
    'pre_bg' : '#272A35',
    'code_bg' : '#272A35',
    'code_text' : '#aebbfe',
}
html_static_path = ['_static']
