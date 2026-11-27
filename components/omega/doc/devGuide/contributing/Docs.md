(omega-dev-docs)=

# Documentation

The Omega documentation is generated using the
[Sphinx](https://www.sphinx-doc.org/en/master/) package and is written in
[MyST](https://myst-parser.readthedocs.io/en/latest/syntax/syntax.html)
format, in `components/omega/doc`. We recommend these
[basic tips](https://myst-parser.readthedocs.io/en/latest/syntax/roles-and-directives.html#roles-directives)
on using MyST in Sphinx.

## Principles

These are the rules a documentation change is reviewed against. The reasons
for them are in the {ref}`documentation design <omega-design-reorg-docs>`.

1. **Organize by the reader's question, not by the code.** A page belongs in
   the guide whose question it answers (see below), whichever source file
   implements it.
2. **One home per fact; everywhere else, link.** Every configuration option,
   equation, class and procedure is described in exactly one place, and
   other pages link to it. Link to Polaris and E3SM documentation rather than
   repeating it.
3. **Docs ship with the code.** A pull request that changes behavior updates
   the affected pages in the same pull request.
4. **Design documents are frozen.** They record the reasoning at the time of
   design and are not updated afterward; the guides describe the code as it
   is.
5. **Generate reference material; write everything else.** Lists that must be
   complete and exact, such as the configuration options, come from one
   machine-readable source.
6. **Math is written to be checked against code.** Technical Guide notation
   stays close to the variable names in the source, and each page links to
   the files that implement it.
7. **Placeholders are explicit, never silent.** A page that has not been
   written says what belongs there and links the issue that tracks it.
8. **Every page has a stable anchor**, prefixed by its guide.
9. **Prefer short, focused pages over long chapters.**
10. **Cite the literature through the bibliography.**

## Where things go

| Guide | Question it answers | Sections |
| ----- | ------------------- | -------- |
| User Guide (`userGuide/`) | How do I set up and run Omega, and what can I control? | `build/`, `run/`, `config/`, `io/`, `model/`, `tutorials/` |
| Developer Guide (`devGuide/`) | How do I work on the code? | `gettingStarted/`, `contributing/`, `framework/`, `model/`, `howto/` |
| Technical Guide (`techGuide/`) | What does the model compute, and why? | `equations/`, `discretization/`, `physics/` |
| Design documents (`design/`) | What were we thinking when we built it? | one directory |

A new feature usually touches all three guides:

- the **User Guide** page for the feature says what a user can set (config
  options, build flags) and what it produces, and links the Technical Guide
  for what it means;
- the **Developer Guide** page describes the classes and how to use them, and
  the unit and Polaris tests that cover the feature;
- the **Technical Guide** page gives the continuous equations and their
  discretized form, with variables close to those in the code.

## Adding a page

Put the page in the directory of the section it belongs to, with a CamelCase
file name, and add it to the toctree in that section's `index.md`. Start the
page with its anchor and a level-1 title:

```markdown
(omega-user-model-frazil)=

# Frazil
```

Anchors start with the guide's prefix: `omega-user-<section>-<page>` in the
User Guide, `omega-dev-<page>` in the Developer Guide, `omega-tech-<page>` in
the Technical Guide and `omega-design-<page>` for design documents. Sphinx
treats anchors without regard to case, so two that differ only in case
collide. Link to other pages with `` {ref}`anchor` `` or `[text](#anchor)`
rather than a relative path to the Markdown file, so that links survive when
pages move.

Images go in the `images/` directory of the guide that uses them.

A page that has not been written yet is a placeholder:

```markdown
(omega-user-quick-start)=

# Quick Start

:::{admonition} Not yet written
:class: warning
This page is a placeholder. It will cover: <one or two sentences>.
Progress is tracked in [E3SM-Project/Omega#NNN](https://github.com/E3SM-Project/Omega/issues/NNN).
:::
```

### Citations

Technical Guide pages cite papers from `techGuide/references.bib`, which the
{ref}`omega-tech-bibliography` page lists. Add an entry there, keyed by the
first author's surname and the year, and cite it with `` {cite:t}`Key` ``
("Ringler et al. (2010)") or `` {cite:p}`Key` `` ("(Ringler et al., 2010)").
Only cited entries appear in the bibliography.

(omega-dev-building-docs)=

## Building the documentation

As long as you have followed the procedure in {ref}`omega-dev-env` for
setting up your development environment, you will already have the packages
available that you need to build the documentation.

From the root of an Omega branch, run the following to build the docs:

```bash
cd components/omega/doc
make clean html
```

You can view the documentation by opening `_build/html/index.html` in your
choice of browsers or copy it to a web portal.

Pull requests are checked with `make html-strict`, which turns warnings into
errors; run it before you push. If you run into errors or warnings related to
changes you have made, check with the development team for help on fixing
them.

### Versioned docs

The published docs keep one directory per version, with a version switcher
in the sidebar. `develop` is published whenever a change to the docs merges
into `develop`, and each release when its `Omega-v*` tag is pushed. To
preview that layout locally, run:

```bash
cd components/omega/doc
make versioned-html
python -m http.server -d _build/html
```

and open <http://localhost:8000/develop/index.html>. The switcher only works
through a web server like this one, since browsers do not let it read the
list of versions from local files. It lists every version built into
`_build/html`; build another with, for example,
`make versioned-html DOCS_VERSION=0.1.0`. `make versioned-html-strict` turns
warnings into errors, and `make clean-versioned-html` removes all versions.

(omega-dev-docs-design)=

## Design documents

A design document records the design of a new capability before it is
implemented. Start from the {doc}`template </design/Template>`, put the new
document in `design/`, and add it to the toctree in `design/index.md`, before
`Template`. Once the design is accepted, the document is not updated to follow
the code; the guides are.
