# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import dataclasses
import os
from datetime import date

import sphinxcontrib.bibtex.plugin
from sphinxcontrib.bibtex.style.referencing import BracketStyle
from sphinxcontrib.bibtex.style.referencing.author_year import (
    AuthorYearReferenceStyle,
)

# -- Project information -----------------------------------------------------

project = "Omega"
copyright = f"{date.today().year}, Energy Exascale Earth System Model Project"
author = "E3SM Development Team"

# The version of the docs: develop or a release version such as 0.1.0. It
# replaces |version| and |release| and names the build in the version
# switcher. The publish workflow sets DOCS_VERSION; see the Makefile.
version = os.environ.get("DOCS_VERSION", "develop")
release = version

master_doc = "index"
language = "en"

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "myst_parser",
    "sphinx_rtd_theme",
    "sphinx.ext.viewcode",
    "sphinx.ext.mathjax",
    "sphinxcontrib.bibtex",
]

numfig = True
# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "README.md",
    "**/README.md",
]

# -- MyST settings ---------------------------------------------------

myst_enable_extensions = ["colon_fence", "deflist", "dollarmath"]
math_number_all = True
math_numfig = False
myst_dmath_allow_labels = True
myst_number_code_blocks = ["typescript"]
myst_heading_anchors = 3
myst_footnote_transition = True
myst_dmath_double_inline = True
myst_enable_checkboxes = True
suppress_warnings = ["myst.header"]
myst_heading_start_level = 1

# -- Bibliography ------------------------------------------------

# Pages cite with {cite:t} ("Ringler et al. (2010)") or {cite:p}
# ("(Ringler et al., 2010)"); techGuide/Bibliography.md lists every cited work.
bibtex_bibfiles = ["techGuide/references.bib"]
bibtex_default_style = "plain"
bibtex_reference_style = "author_year_round"


def _round_brackets():
    return BracketStyle(left="(", right=")", sep="; ")


@dataclasses.dataclass
class _AuthorYearRoundStyle(AuthorYearReferenceStyle):
    """Author-year citations in round brackets, separated by semicolons"""

    bracket_parenthetical: BracketStyle = dataclasses.field(
        default_factory=_round_brackets
    )
    bracket_textual: BracketStyle = dataclasses.field(
        default_factory=_round_brackets
    )
    bracket_author: BracketStyle = dataclasses.field(
        default_factory=_round_brackets
    )
    bracket_label: BracketStyle = dataclasses.field(
        default_factory=_round_brackets
    )
    bracket_year: BracketStyle = dataclasses.field(
        default_factory=_round_brackets
    )


sphinxcontrib.bibtex.plugin.register_plugin(
    "sphinxcontrib.bibtex.style.referencing",
    "author_year_round",
    _AuthorYearRoundStyle,
)

# -- HTML output -------------------------------------------------

html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

# Used by _templates/layout.html to select this version in the switcher
html_context = {"current_version": version}
