# coding: utf-8
# Copyright (c) 2018-2019, Taku MURAKAMI. All rights reserved.
# Distributed under the terms of the BSD 3-clause License.

from setuptools import setup
from setuptools import find_packages


with open("README.md") as file:
    readme = file.read()

with open('LICENSE') as file:
    license = file.read()

setup(
    name="pythroughput",
    version="1.0",
    description="Python module to perform high-throughput first-principles calculation in 'Xenonpy' package.",
    long_description=readme,
    author="Taku MURAKAMI",
    author_email="murakami.taku.17@shizuoka.ac.jp",
    url="https://github.com/murakami17/pythroughput",
    license=license,
    packages=find_packages(exclude=("tests", "docs"))
)
