#!/usr/bin/env python

import sys
import os

from setuptools import setup, Extension, find_packages

# from distutils.ccompiler import CCompiler
# from distutils.unixccompiler import UnixCCompiler
# from distutils.msvccompiler import MSVCCompiler

# local_path = os.path.dirname(os.path.abspath(__file__))

local_path = os.path.dirname(os.path.abspath(__file__))
# change compiler
# os.environ["CC"] = "clang++"

# os.environ["CC"] = "g++"
# os.environ["CC"] = "cl"

# turn off warnings raised by Minuit and generated Cython code that need
# to be fixed in the original code bases of Minuit and Cython
# compiler_opts = {
#     CCompiler: {},
#     UnixCCompiler: {
#         "extra_compile_args": [
#             "-std=c++11",
#             "-O3"
#         ]
#     },
#     MSVCCompiler: {"extra_compile_args": ["/EHsc", "/O2"]},
# }

# C/C++ extension of gtracr module
# contains magnetic field and trajactory evaluation as core of the code
_libgtracr = Extension(
    '_libgtracr',
    sources=[
        "gtracr/lib/src/TrajectoryTracer.cpp",
        "gtracr/lib/src/uTrajectoryTracer.cpp", "gtracr/lib/src/igrf.cpp",
        "gtracr/lib/src/pybind11_wrapper.cpp"
    ],
    language='c++',
    include_dirs=['gtracr/lib/include']
)
# extra_compile_args=[
#     "-fsave-optimization-record=yaml", "-fdiagnostics-show-hotness",
#     "-Rpass=.*", "-foptimization-record-file=" +
#     os.path.join(local_path, "opt_reports", "_gtracr.opt.yaml")
# ]
# extra_compile_args=["-O3"]),
# extra_compile_args=["/O2"]),


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

exclude_dirs = ["*.tests", "*.tests.*", "tests.*", "tests",
                "*.data", "*.data.*", "data.*", "data"]


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
    # packages=['gtracr', 'gtracr.scripts', 'gtracr.tests',
    #           'gtracr.lib', 'gtracr.lib.magnetic_field'],

    packages=find_packages(
        exclude=exclude_dirs),  # finds dir with __init__.py, excludes tests
    incluse_package_data=True,
    ext_modules=[_libgtracr],
    install_requires=[
        'scipy',
        'numpy' 
    ],
    extras_require=extras_require,
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3', 'Topic :: Software Development',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Physics',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Development Status :: 1 - Planning', 'Natural Language :: English',
        'License :: OSI Approved :: BSD License'
    ],
    python_requires='>=3.0',
)
