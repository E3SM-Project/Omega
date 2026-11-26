(omega-dev-building-docs)=

# Building the Documentation

As long as you have followed the procedure in {ref}`omega-dev-env` for
setting up your development environment, you will already have the packages
available that you need to build the documentation.

From the root of an Omega branch, run the following script to build the docs:

```bash
cd components/omega/doc
make clean html
```

You can view the documentation by opening `_build/html/index.html` in your
choice of browsers or copy it to a web portal.

If you run into errors or warnings related to changes you have made when you
build the docs, check with the development team for help on fixing them.

## Versioned docs

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
