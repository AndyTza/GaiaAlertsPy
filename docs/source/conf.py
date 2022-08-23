# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# -- Project information -----------------------------------------------------

project = "GaiaAlertPy"
copyright = "2022, GaiaAlertsPy"
author = "Anastasios Tzanidakis"

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx_copybutton",
    "nbsphinx",
    "sphinx_automodapi.automodapi",
]
numpydoc_show_class_members = False

# Add any paths that contain templates here, relative to this directory.

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# hide input prompts in notebooks
nbsphinx_prompt_width = "0"

#Theme
html_theme = "furo"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]


# hide coppy button on outputs
copybutton_selector = "div:not(.output_area) > div.highlight pre"