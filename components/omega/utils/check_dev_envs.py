#!/usr/bin/env python3
"""
Check that dev-conda.txt and pixi.toml list the same packages.

pixi.toml is the preferred development environment and dev-conda.txt is a
fallback for developers who use conda. This script compares the packages in
pixi's default environment with those in dev-conda.txt, including their
version constraints. Python is skipped because only pixi.toml pins it.

Usage (from components/omega):
$ utils/check_dev_envs.py
"""

import argparse
import re
import sys
from pathlib import Path

import tomllib

SKIP = {"python"}


def main():
    """
    Compare dev-conda.txt and pixi.toml and exit with an error if they differ
    """
    omega_dir = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument(
        "--conda", type=Path, default=omega_dir / "dev-conda.txt",
        help="the conda spec file",
    )  # fmt: skip
    parser.add_argument(
        "--pixi", type=Path, default=omega_dir / "pixi.toml",
        help="the pixi manifest",
    )  # fmt: skip
    args = parser.parse_args()

    conda = read_conda_spec(args.conda)
    pixi = read_pixi_default_env(args.pixi)

    problems = []
    for name in sorted(conda.keys() - pixi.keys()):
        problems.append(f"{name} is in {args.conda.name} but not pixi.toml")
    for name in sorted(pixi.keys() - conda.keys()):
        problems.append(f"{name} is in pixi.toml but not {args.conda.name}")
    for name in sorted(conda.keys() & pixi.keys()):
        if conda[name] != pixi[name]:
            problems.append(
                f"{name} is '{conda[name] or '*'}' in {args.conda.name} but "
                f"'{pixi[name] or '*'}' in pixi.toml"
            )

    if problems:
        print(f"{args.conda.name} and pixi.toml differ:")
        for problem in problems:
            print(f"  {problem}")
        sys.exit(1)

    print(f"{args.conda.name} and pixi.toml list the same packages")


def read_conda_spec(path):
    """
    Read package names and version constraints from a conda spec file

    Parameters
    ----------
    path : pathlib.Path
        The conda spec file

    Returns
    -------
    packages : dict
        Version constraints (empty for none) keyed by package name
    """
    packages: dict[str, str] = {}
    for line in path.read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        match = re.match(r"([A-Za-z0-9_.\-]+)(.*)", line)
        if match is None:
            msg = f"Could not parse '{line}' in {path}"
            raise ValueError(msg)
        name, constraint = match.groups()
        _add_package(packages, name, constraint, path)
    return packages


def read_pixi_default_env(path):
    """
    Read package names and version constraints of pixi's default environment

    Parameters
    ----------
    path : pathlib.Path
        The pixi manifest

    Returns
    -------
    packages : dict
        Version constraints (empty for none) keyed by package name
    """
    with open(path, "rb") as f:
        manifest = tomllib.load(f)

    tables = [manifest.get("dependencies", {})]
    default_env = manifest.get("environments", {}).get("default", [])
    if isinstance(default_env, dict):
        default_env = default_env.get("features", [])
    for feature in default_env:
        feature_table = manifest.get("feature", {}).get(feature, {})
        tables.append(feature_table.get("dependencies", {}))

    packages: dict[str, str] = {}
    for table in tables:
        for name, spec in table.items():
            if isinstance(spec, dict):
                spec = spec.get("version", "*")
            _add_package(packages, name, spec, path)
    return packages


def _add_package(packages, name, constraint, path):
    """
    Add a package with a normalized name and constraint, skipping python
    """
    name = name.lower().replace("_", "-")
    if name in SKIP:
        return
    constraint = re.sub(r"\s+", "", constraint)
    if constraint == "*":
        constraint = ""
    if name in packages:
        msg = f"{name} is listed more than once in {path}"
        raise ValueError(msg)
    packages[name] = constraint


if __name__ == "__main__":
    main()
