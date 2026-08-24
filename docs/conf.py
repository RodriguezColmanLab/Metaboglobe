# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

from datetime import date
current_year = date.today().year
copyright_years = f'2026-{current_year}' if current_year > 2026 else f'{current_year}'  # Can be simplified in 2027
project = 'Metaboglobe'
copyright = f'{copyright_years}, Rodríguez Colman lab, UMC Utrecht.'
author = 'Rutger Kok'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'myst_parser',  # Markdown support
]


templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'pydata_sphinx_theme'
html_static_path = ['_static']
html_theme_options = {
  "external_links": [
      {"name": "Source Code", "url": "https://github.com/RodriguezColmanLab/Metaboglobe"}
  ]
}

# -- Options for Markdown input -------------------------------------------------
myst_enable_extensions = [
    "colon_fence",
]
