#!/usr/bin/env python

import sys
import os

from setuptools import setup, Extension, find_packages

local_path = os.path.dirname(os.path.abspath(__file__))


setup(
    name='GSimTraCR',
    version='1.0.0',
    description=
    'A GPU-based simulation that tracks cosmic rays using real geomagnetic field data from IGRF. This project is purely Pythonic (for now).',
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
)

# extra_compile_args=[    # calls error when run by standard C++ compiler (Windows)
#     '-std=c++14',   # require C++14 or higher
#     '-Wno-unused-function',  # do not raise error even when function is unused
#     '-Wno-write-strings',
# ],