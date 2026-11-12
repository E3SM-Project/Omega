(omega-dev-env)=

# Development Environment

As a developer, if you want to lint the code, build the documentation locally
or run the `omega_buildnml` unit tests, you will need your own environment with
the necessary packages.

We recommend [pixi](https://pixi.sh), which creates the environment from
`components/omega/pixi.toml`. Omega's CI also uses this environment. If you
cannot use pixi, {ref}`omega-dev-conda-env` describes how to create a conda
environment from `components/omega/dev-conda.txt` instead.

(omega-dev-install-pixi)=

## Installing pixi

If you have not already done so, install pixi:
```bash
curl -fsSL https://pixi.sh/install.sh | sh
```
This installs `pixi` in `${HOME}/.pixi/bin` and adds that directory to your
`PATH` in your `.bashrc`. Start a new shell (or `source ~/.bashrc`) before
using it. Pixi is a single small executable, so it can stay in your home
directory even where your home quota is small.

(omega-dev-pixi-cache)=

## Keeping pixi off your home quota

Before you create any pixi environments on an HPC machine, decide where pixi
will keep its packages.

Pixi downloads each package once into a package cache, and then hardlinks the
files from the cache into each environment. Hardlinks only work when the cache
and the environment are on the same file system. Otherwise, pixi copies the
packages into every environment. The Omega environment is about 2 GB, and
every branch or worktree gets its own environment, so copies add up quickly.
Home directories on HPC machines usually have small quotas.

To avoid this, add a line like the following to your `.bashrc` on each
machine, pointing at a directory with plenty of space:
```bash
export PIXI_CACHE_DIR=/path/with/space/pixi/cache
```
Set `PIXI_CACHE_DIR` explicitly even if the default location would have room.
On some network file systems, pixi otherwise caches packages inside each
workspace, so every branch or worktree gets its own copy of every package.

Pixi puts an environment in `components/omega/.pixi` of the branch you create
it from. If your Omega branches are not on the same file system as
`PIXI_CACHE_DIR` (for example, if your branches are in your home directory but
the cache is on a project file system), also tell pixi to put its environments
next to the cache:
```bash
pixi config set --global detached-environments /path/with/space/pixi/envs
```
You only need to run this once on each machine. It is saved in
`${HOME}/.pixi/config.toml`.

For example, on Perlmutter, the home directory is too small for Omega branches,
so we recommend cloning Omega somewhere under `/global/cfs/cdirs/e3sm/${USER}`.
Then put the cache on the same file system:
```bash
# in ~/.bashrc
export PIXI_CACHE_DIR=/global/cfs/cdirs/e3sm/${USER}/pixi/cache
```
If you do keep Omega branches in your home directory, this cache alone would
make pixi copy packages into each environment in your home directory. In that
case, also put the environments on CFS:
```bash
# once
pixi config set --global detached-environments /global/cfs/cdirs/e3sm/${USER}/pixi/envs
```

Avoid scratch file systems that are purged. If files in an environment are
purged, the environment stops working and you will need to recreate it.

On a laptop or workstation, the default locations are usually fine.

(omega-dev-create-env)=

## Creating the environment

From the root of an Omega branch, run:
```bash
cd components/omega/
pixi install
```
Pixi solves for the packages in `pixi.toml` and installs them. It also writes a
`pixi.lock` file recording exactly what it installed. Omega does not commit
this file to the repository (it is in `.gitignore`), so each branch gets the
latest compatible packages when its environment is created.

(omega-dev-activate-env)=

## Activating the environment

To activate the environment, run:
```bash
cd components/omega/
pixi shell
```
This starts a new shell with the environment activated. You can then move
around the repository (or anywhere else) and use the tools in the environment.
Type `exit` to leave the shell. You can also run `pixi shell` from the root of
the branch with:
```bash
pixi shell --manifest-path components/omega/pixi.toml
```
If you have not run `pixi install` yet, `pixi shell` will create the
environment first.

The first time you set up the environment, you need to initiate `pre-commit`:
```bash
pre-commit install
```
This needs to be done once per environment, so re-run it as necessary if you
re-create the development environment.

Please activate this environment each time you commit code. This will ensure
that `pre-commit` and the associated linting utilities are available and that
they check your code as it is committed (rather than requiring fix-up commits
later on).

(omega-dev-update-env)=

## Updating the environment

If `pixi.toml` changes (for example, after you merge or rebase onto a newer
`develop`), pixi updates the environment the next time you run `pixi shell` or
`pixi install`.

To update the packages to their latest compatible versions without any change
to `pixi.toml`, run:
```bash
cd components/omega/
pixi update
```

(omega-dev-change-env)=

## Changing the packages

If you add, remove or change the version constraint of a package in
`pixi.toml`, make the same change in `dev-conda.txt`. CI runs
`components/omega/utils/check_dev_envs.py` to make sure the two files list the
same packages and version constraints. You can run it yourself before you
commit:
```bash
cd components/omega/
utils/check_dev_envs.py
```

(omega-dev-conda-env)=

## Fallback: conda environment

If you cannot use pixi, you can create a
[conda](https://conda.io/projects/conda/en/latest/index.html) environment with
the same packages instead. Unlike the pixi environment, this environment is not
built in CI, though CI does check that `dev-conda.txt` lists the same packages
as `pixi.toml`.

(omega-dev-install-miniforge3)=

### Installing Miniforge3

If you have not already done so, you will install
[Miniforge3](https://github.com/conda-forge/miniforge#miniforge3):
```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
```
If you have the space to do so, we recommend installing Miniforge3 somewhere in
our home directory.  Typically, this file system will be faster than on a
scratch space spaces.  Also, a conda environment will be rendered useless if
its contents get purged.  However, conda environments each require several GB
of space, which may quickly use up your disk quota so you will need to weigh
these trade-offs for yourself.

In what follows, we assume you have Miniforge3 installed in
`${HOME}/miniforge3`.

You may wish to skip the step in the Miniforge3 installation where it
wants to modify your `.bashrc`.  If you allow it to, it will add activation of
the `conda` command and the `base` conda environment, and this may not be
something you want.

If you opt not to allow Miniforge3 to update your `.bashrc`, it may be helpful
to add an alias like the following to your `.bashrc`:

```bash
alias conda_base='
source ${HOME}/miniforge3/etc/profile.d/conda.sh
conda activate'
```

Anytime you need `conda` command, you can first run `conda_base`to get it.

(omega-dev-create-conda-env)=

### Creating the Conda Environment

With Miniforge3 installed and the `conda` command available, you can create the
Omega development environment, starting from the root
of an Omega branch, with:

```bash
cd components/omega/
conda create -n omega_dev --file dev-conda.txt
conda activate omega_dev
pre-commit install
```

(You can reuse `omega_dev` for other branches as long as `dev-conda.txt` has
not changed between branches.)

(omega-dev-update-conda-env)=

### Updating the Conda Environment

If you need to update the `omega_dev` conda environment, the easiest way to
do that is just to recreate it by following the same instructions as above.
You will be prompted about whether you want to delete the old environment
unless you include the `-y` flag (to say "yes" to all prompts).
