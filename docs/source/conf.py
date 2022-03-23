# Configuration file for the Sphinx documentation builder.
# using commands from https://gist.github.com/matteoferla/ba72ab12a9e5f690277e2e88551773aa
# modified for readthedocs
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# ``.readthedocs.yaml` installs it.
# import os
# import sys
# sys.path.insert(0, os.path.abspath("../../"))

# -- Project information -----------------------------------------------------

project = 'smallworld_api'
copyright = '2022, Matteo Ferla'
author = 'Matteo Ferla'
github_username = 'matteoferla'
github_repository = 'Python_SmallWorld_API'


# -- General configuration ---------------------------------------------------

extensions = [
    'readthedocs_ext.readthedocs',
    'sphinx.ext.viewcode',
    'sphinx.ext.todo',
    #'sphinx_toolbox.more_autodoc',
    'sphinx.ext.autodoc',
    #'sphinx.ext.imgconverter',
    #'m2r'
]

html_static_path = []
templates_path = ['_templates']
#source_suffix = ['.rst', '.md']
always_document_param_types = True
typehints_defaults = 'braces'

from m2r2 import parse_from_file  # noqa

for markdown_filename, srt_filename in {'../../README.md': 'readme.rst',
                                        '../../export_note.md': 'export_note.rst'}.items():
    with open(srt_filename, 'w') as fh:
        fh.write(parse_from_file(markdown_filename))

language = 'en'
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
# html_theme = 'alabaster'
todo_include_todos = True

def skip(app, what, name, obj, would_skip, options):
    if name in ( '__init__',):
        return False
    return would_skip

def setup(app):
    app.connect('autodoc-skip-member', skip)