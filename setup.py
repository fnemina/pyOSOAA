#!/usr/bin/env python
# coding=utf-8
# This file is part of pyOSOAA.
#
# Copyright 2019 Francisco Nemiña

import os
from setuptools import setup

PROJECT_ROOT = os.path.dirname(__file__)


def read_file(filepath, root=PROJECT_ROOT):
    """
    Return the contents of the specified `filepath`.
    * `root` is the base path and it defaults to the `PROJECT_ROOT` directory.
    * `filepath` should be a relative path, starting from `root`.
    """
    with open(os.path.join(root, filepath), encoding="utf-8") as fd:
        text = fd.read()
    return text


LONG_DESCRIPTION = read_file("README.md")
SHORT_DESCRIPTION = (
    "pyOSOAA is a python interface for the Ocean Successive "
    "Orders with Atmosphere - Advanced (OSOAA) radiative transfer."
)

setup(
    name="pyOSOAA",
    version="1.5.0",
    author="Francisco Nemiña",
    author_email="fnemina@conae.gov.ar",
    description=SHORT_DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    url="https://github.com/fnemina/pyOSOAA",
    project_urls={
        "Bug Tracker": "https://github.com/fnemina/pyOSOAA/issues",
    },
    packages=["pyOSOAA"],
    python_requires=">=3.8",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Atmospheric Science",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
)
