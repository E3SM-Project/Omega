#!/usr/bin/env python3
"""
Write shared/versions.json, the list of documentation versions for the
version switcher.

Each version of the docs is built into its own directory under a base
directory, next to a shared/ directory that holds the switcher script:

    <base>/develop/
    <base>/0.1.0/
    <base>/shared/version-switcher.js
    <base>/shared/versions.json

A version directory is any directory in <base> that is the root of a Sphinx
HTML build. develop comes first, then release versions, newest first, then any
other versions alphabetically. URLs are relative to versions.json.

Usage (from components/omega/doc):
$ ./generate_versions_json.py --base-dir _build/html
"""

import argparse
import json
import re
from pathlib import Path

DEFAULT_VERSION = "develop"

# Semantic versions such as 1.2.3 and 0.1.0-alpha.1
RELEASE_RE = re.compile(r"^(\d+)\.(\d+)\.(\d+)(?:-([0-9A-Za-z.-]+))?$")


def main():
    """
    Find the version directories and write their list to versions.json
    """
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument(
        "--base-dir", type=Path, required=True,
        help="the directory holding one directory per version",
    )  # fmt: skip
    args = parser.parse_args()

    base_dir = args.base_dir
    if not base_dir.is_dir():
        msg = f"Base directory {base_dir} does not exist"
        raise FileNotFoundError(msg)

    versions = sort_versions(
        path.name for path in base_dir.iterdir() if is_version_dir(path)
    )
    if not versions:
        msg = f"No documentation versions found in {base_dir}"
        raise ValueError(msg)

    entries = [
        {"version": version, "url": f"../{version}/index.html"}
        for version in versions
    ]

    shared_dir = base_dir / "shared"
    shared_dir.mkdir(exist_ok=True)
    with open(shared_dir / "versions.json", "w") as f:
        json.dump(entries, f, indent=2)
        f.write("\n")

    print(f"Wrote {shared_dir / 'versions.json'}: {', '.join(versions)}")


def is_version_dir(path):
    """
    Whether a path is the root of a Sphinx HTML build

    Parameters
    ----------
    path : pathlib.Path
        A path in the base directory

    Returns
    -------
    is_version : bool
        True if the path is a directory with a Sphinx inventory and index
    """
    return (
        path.is_dir()
        and (path / "objects.inv").is_file()
        and (path / "index.html").is_file()
    )


def sort_versions(versions):
    """
    Sort versions for the version switcher

    Parameters
    ----------
    versions : iterable of str
        The version names

    Returns
    -------
    sorted_versions : list of str
        develop, then releases from newest to oldest, then other versions
        alphabetically
    """
    versions = set(versions)
    releases = sorted(
        (v for v in versions if RELEASE_RE.match(v)),
        key=release_key,
        reverse=True,
    )
    others = sorted(
        v for v in versions if v != DEFAULT_VERSION and v not in releases
    )
    default = [DEFAULT_VERSION] if DEFAULT_VERSION in versions else []
    return default + releases + others


def release_key(version):
    """
    A sort key that orders release versions by semantic-version precedence

    A pre-release such as 1.0.0-alpha.1 comes before the release 1.0.0.
    Pre-release identifiers compare numerically if they are numbers and
    alphabetically otherwise, so alpha.2 < alpha.10 < beta.1 < rc.1.

    Parameters
    ----------
    version : str
        A release version matching ``RELEASE_RE``

    Returns
    -------
    key : tuple
        The sort key
    """
    match = RELEASE_RE.match(version)
    if match is None:
        msg = f"{version} is not a release version"
        raise ValueError(msg)
    major, minor, patch, pre = match.groups()
    numbers = (int(major), int(minor), int(patch))
    if pre is None:
        return numbers + (1, ())
    identifiers = tuple(
        (0, int(part), "") if part.isdigit() else (1, 0, part)
        for part in pre.split(".")
    )
    return numbers + (0, identifiers)


if __name__ == "__main__":
    main()
