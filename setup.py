#!/usr/bin/env python

import sys
import os

from setuptools import setup, Extension, find_packages

# local_path = os.path.dirname(os.path.abspath(__file__))


with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name='gtracr',
    version='1.0.0',
    description=
    'A GPU-based simulation that tracks cosmic rays.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/kwat0308/gtracr",
    packages=find_packages(exclude='tests'),
    classifiers=[
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Topic :: Software Development',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Physics',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Development Status :: 1 - Planning',
        'Environment :: Console',
        'Natural Language :: English',
    ],
    python_requires='>=3.0',
)

# extra_compile_args=[    # calls error when run by standard C++ compiler (Windows)
#     '-std=c++14',   # require C++14 or higher
#     '-Wno-unused-function',  # do not raise error even when function is unused
#     '-Wno-write-strings',
# ],