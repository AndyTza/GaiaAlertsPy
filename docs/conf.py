# conf.py for GaiaAlertsPy documentation
import os
import sys
sys.path.insert(0, os.path.abspath('.'))

# Project information
project = 'GaiaAlertsPy'
author = 'Andy Tzanidakis'
release = '0.1'

# General configuration
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'nbsphinx'
]

templates_path = ['_templates']
exclude_patterns = []

# HTML output
html_theme = 'pydata_sphinx_theme'

html_static_path = ['_static']

html_css_files = [
    'css/custom.css',
]

html_context = {"default_mode": "light"}

html_theme_options = {
    'canonical_url': '',
    'default_mode': 'light',
    'display_version': True,
    'navbar_end': ['navbar-icon-links.html'],
    'logo': {
        'text': project,
    },
    'light_mode': True,  # Enable light mode by default
    'dark_mode': False,   # Enable dark mode option
}

