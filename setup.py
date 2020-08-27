#!/usr/bin/env python

import sys
import os

from setuptools import setup, Extension, find_packages

local_path = os.path.dirname(os.path.abspath(__file__))

# C/C++ extension of gtracr module
# contains magnetic field and trajactory evaluation as core of the code
libgtracr = Extension(
    'gtracr.lib.libgtracr',
    sources=[
        "gtracr/lib/src/TrajectoryTracer.cpp",
        "gtracr/lib/src/uTrajectoryTracer.cpp", "gtracr/lib/src/igrf.cpp",
        "gtracr/lib/src/pybind11_wrapper.cpp"
    ],
    language='c++',
    include_dirs=['gtracr/lib/include']
)


# This method is adopted from MCEq https://github.com/afedynitch/MCEq
# which is adopted from iMinuit https://github.com/scikit-hep/iminuit
# Getting the version number at this point is a bit tricky in Python:
# https://packaging.python.org/en/latest/development.html#single-sourcing-the-version-across-setup-py-and-your-project
# This is one of the recommended methods that works in Python 2 and 3:
def get_version():
    version = {}
    with open("gtracr/version.py") as fp:
        exec(fp.read(), version)
    return version['__version__']


__version__ = get_version()

extras_require = {
    "test": ["pytest", "pytest-benchmark"],
    "examples": ["matplotlib", "mpld3"]
}

# exclude_dirs = ["*.tests", "*.tests.*", "tests.*", "tests",
#                 "*.data", "*.data.*", "data.*", "data"]


# open README for the long descreption of the code
with open("README.md", "r") as f:
    long_description = f.read()

# main setup configuration
setup(
    name='gtracr',
    version=__version__,
    description='A GPU-based simulation that tracks cosmic rays from any location on Earth.',
    author='Keito Watanabe',
    author_email='k.wat8973@gmail.com',
    long_description=long_description,
    long_description_content_type="text/markdown",
    license='BSD 3-Clause License',
    url="https://github.com/kwat0308/gtracr",
    packages=['gtracr', 'gtracr.scripts', 'gtracr.tests',
              'gtracr.lib'],
    # include_package_data=True,
    package_data={
        'gtracr': ['data/**'],
    },
    ext_modules=[libgtracr],
    install_requires=[
        'scipy',
        'numpy',
        'datetime'
    ],
    extras_require=extras_require,
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3', 'Topic :: Software Development',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Physics',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Development Status :: 4 - Beta', 'Natural Language :: English',
        'License :: OSI Approved :: BSD License'
    ],
    python_requires='>=3.0',
)
