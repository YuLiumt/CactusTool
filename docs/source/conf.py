# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('../..'))

package_path = os.path.abspath('../..')
os.environ['PYTHONPATH'] = ':'.join((package_path, os.environ.get('PYTHONPATH', '')))

# autodoc_mock_imports = ['numpy', 're', 'os', 'h5py', 'bz2', 'gzip']
import importlib
autodoc_mock_imports = []
for mod in ['numpy', 're', 'os', 'h5py', 'bz2', 'gzip']:
    try:
        importlib.import_module(mod)
    except ImportError:
        autodoc_mock_imports.append(mod)


import sphinx_rtd_theme
# -- Project information -----------------------------------------------------

project = 'CactusTool'
copyright = '2020, Yu Liu'
author = 'Yu Liu'

# The full version, including alpha/beta/rc tags
release = '0.0.1'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.githubpages',
    'sphinx.ext.todo',
    'sphinx.ext.viewcode',
    'sphinx.ext.graphviz',
    'sphinx.ext.mathjax',
    "sphinx_rtd_theme",
    'sphinx.ext.autosummary',
    'nbsphinx',
    'IPython.sphinxext.ipython_console_highlighting',
    'recommonmark',
    # 'jupyter_sphinx.execute',
]

nbsphinx_allow_errors = True
# nbsphinx_execute = 'always'

autosummary_generate = True

todo_include_todos = True
# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = [
    '_build',
    '**.ipynb_checkpoints',
]

language = 'en'

master_doc = 'index'
source_suffix = '.rst'
# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = "sphinx_rtd_theme"
html_theme = "nature"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']