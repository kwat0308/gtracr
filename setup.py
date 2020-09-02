#!/usr/bin/env python

import sys
import os
import platform
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from distutils.ccompiler import CCompiler
from distutils.unixccompiler import UnixCCompiler
from distutils.msvccompiler import MSVCCompiler
import distutils.ccompiler

local_path = os.path.dirname(os.path.abspath(__file__))

# C/C++ extension of gtracr module
# contains magnetic field and trajactory evaluation as core of the code
libgtracr = Extension('gtracr.lib._libgtracr',
                      sources=[
                          "gtracr/lib/src/TrajectoryTracer.cpp",
                          "gtracr/lib/src/uTrajectoryTracer.cpp",
                          "gtracr/lib/src/igrf.cpp",
                          "gtracr/lib/src/pybind11_wrapper.cpp"
                      ],
                      language='c++',
                      include_dirs=['gtracr/lib/include', 'gtracr/lib/extern'])

'''
The below settings were obtained from the Iminuit package from scikit-HEP:
https://github.com/scikit-hep/iminuit 
'''


extra_flags = []
if bool(os.environ.get("COVERAGE", False)):
    extra_flags += ["--coverage"]
if platform.system() == "Darwin":
    extra_flags += ["-stdlib=libc++"]

# turn off warnings raised by Minuit and generated Cython code that need
# to be fixed in the original code bases of Minuit and Cython
compiler_opts = {
    CCompiler: {},
    UnixCCompiler: {
        "extra_compile_args": [
            "-std=c++11",
            "-Wno-shorten-64-to-32",
            "-Wno-parentheses",
            "-Wno-unused-variable",
            "-Wno-sign-compare",
            "-Wno-cpp",  # suppresses #warnings from numpy
            "-Wno-deprecated-declarations",
        ]
        + extra_flags,
        "extra_link_args": extra_flags,
    },
    MSVCCompiler: {"extra_compile_args": ["/EHsc"]},
}


class SmartBuildExt(build_ext):
    def build_extensions(self):
        c = self.compiler
        opts = [v for k, v in compiler_opts.items() if isinstance(c, k)]
        for e in self.extensions:
            for o in opts:
                for attrib, value in o.items():
                    getattr(e, attrib).extend(value)

        build_ext.build_extensions(self)

# Getting the version number at this point is a bit tricky in Python:
# https://packaging.python.org/en/latest/development.html#single-sourcing-the-version-across-setup-py-and-your-project
# This is one of the recommended methods that works in Python 2 and 3:
def get_version():
    version = {}
    with open("gtracr/version.py") as fp:
        exec(fp.read(), version)
    return version['__version__']


__version__ = get_version()

extras_require = {"examples": ["matplotlib", "plotly"]}

# exclude_dirs = ["*.tests", "*.tests.*", "tests.*", "tests",
#                 "*.data", "*.data.*", "data.*", "data"]

# open README for the long descreption of the code
with open("README.md", "r") as f:
    long_description = f.read()

# main setup configuration
setup(
    name='gtracr',
    version=__version__,
    description=
    'A GPU-based simulation that tracks cosmic rays from any location on Earth.',
    author='Keito Watanabe',
    author_email='k.wat8973@gmail.com',
    long_description=long_description,
    long_description_content_type="text/markdown",
    license='BSD 3-Clause License',
    url="https://github.com/kwat0308/gtracr",
    packages=['gtracr', 'gtracr.scripts', 'gtracr.tests', 'gtracr.lib'],
    # include_package_data=True,
    package_data={
        'gtracr': ['data/**'],
    },
    ext_modules=[libgtracr],
    cmdclass={"build_ext": SmartBuildExt},
    install_requires=[
        'scipy',
        'numpy',
        'datetime',
        'tqdm'
    ],
    extras_require=extras_require,
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3', 'Topic :: Software Development',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Physics',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers', 'Development Status :: 4 - Beta',
        'Natural Language :: English', 'License :: OSI Approved :: BSD License'
    ],
    python_requires='>=3.0'
)
