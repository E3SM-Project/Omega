(omega-dev-linting)=

# Linting Code

## Enabling Lint Checks

First, please follow the procedure in {ref}`omega-dev-conda-env` to set up a
conda environment for linting the code and building the docs.

Once on each machine you are developing on, before committing code, please run:
```bash
pre-commit install
```

After this, your code will be linted with all the tools discussed below as part
of each call to `git commit`. Any modified files will also be linted when CI
(GitHub Actions) runs on each pull request and merge to the `develop` branch.

You can lint the full code using all the tools with:
```bash
pre-commit run --all-files
``````

## Bypassing Linting

If you wish to commit without first checking the code for lint (e.g. you plan
to fix the code in a later commit), run:
```bash
git commit --no-verify ...
```
with the same arguments you would normally use for `git commit`.


## Linting C++ Code

C++ code is formatted with
[clang-format](https://clang.llvm.org/docs/ClangFormat.html) using the style in
`components/omega/.clang-format`.

The `pre-commit` hook installs its own copy of `clang-format` at the version
pinned in `.pre-commit-config.yaml`, so the result does not depend on which
`clang-format` is on your path. Different major versions of `clang-format`
format some code differently. If your editor formats on save, point it at the
`clang-format` in the `omega_dev` environment, which is pinned to the same
LLVM release (in VS Code, for example, set `C_Cpp.clang_format_path`).

You can run the formatter on its own with:
```bash
pre-commit run clang-format --all-files
```
You can specify one or more files instead of `--all-files`.

When the pinned version changes, the code is reformatted in a single commit
that is listed in `.git-blame-ignore-revs` at the top of the repository. To
have `git blame` skip such commits, run this once in your clone:
```bash
git config blame.ignoreRevsFile .git-blame-ignore-revs
```

## Linting Python Code

The tools used to lint python code include
[flynt](https://github.com/ikamensh/flynt) for enforcing string formatting
with f-strings, [ruff](https://docs.astral.sh/ruff/) for sorting imports and
enforcing the [PEP8](https://peps.python.org/pep-0008/) style guide for python
(`ruff-check`) as well as for formatting the code (`ruff-format`), and
[mypy](https://mypy-lang.org/) for performing variable type checking.

`ruff` is configured in `components/omega/ruff.toml` and `mypy` in
`components/omega/mypy.cfg`.

You can run these tools individually if you need to:
```bash
pre-commit run flynt --all-files
pre-commit run ruff-check --all-files
pre-commit run ruff-format --all-files
pre-commit run mypy --all-files
```
You can specify one more more files instead of `--all-files`.

## Linting Fortran Code

The tool used to lint Fortran code is
[fortitude](https://fortitude.readthedocs.io/en/stable/).

You can run  this tool individually if you need to:
```
pre-commit run fortitude --all-files
```
You can specify one more more files instead of `--all-files`.

## Updating the linting pacakge

To update the linting packages, you just need to recreate the development
conda environment.  See {ref}`omega-dev-update-conda-env`.
