#!/usr/bin/env python

import sys
import os

from setuptools import setup, Extension, find_packages

# local_path = os.path.dirname(os.path.abspath(__file__))

local_path = os.path.dirname(os.path.abspath(__file__))
# change compiler
os.environ["CC"] = "clang++"

# os.environ["CC"] = "g++"
# os.environ["CC"] = "cl"


# obtained from: https://github.com/pybind/python_example/blob/master/setup.py
class get_pybind_include(object):
    """Helper class to determine the pybind11 include path
    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """
    def __str__(self):
        import pybind11
        return pybind11.get_include()


ext_module = [
    Extension(
        '_gtracr',
        sources=[
            "gtracr/lib/src/TrajectoryTracer.cc",
            "gtracr/lib/wrapper/TrajectoryTracerWrapper.cc",
            "gtracr/lib/src/Location.cc",
            "gtracr/lib/wrapper/LocationWrapper.cc",
            "gtracr/lib/src/Particle.cc",
            "gtracr/lib/wrapper/ParticleWrapper.cc",
            # "gtracr/lib/src/MagneticField.cc"
        ],
        #   library_dirs=["gtracr/lib"],
        language='c++',
        include_dirs=[get_pybind_include(), 'gtracr/lib/include'],
        extra_compile_args=[
            "-fsave-optimization-record=yaml", "-fdiagnostics-show-hotness",
            "-Rpass=.*", "-Weverything"
        ]),
]

with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name='gtracr',
    version='0.1.0',
    description=
    'A GPU-based simulation that tracks cosmic rays from any location on Earth.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/kwat0308/gtracr",
    packages=find_packages(exclude='tests'),
    ext_modules=ext_module,
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