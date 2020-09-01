#!/usr/bin/env python

import sys
import os

from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
import setuptools

local_path = os.path.dirname(os.path.abspath(__file__))

# C/C++ extension of gtracr module
# contains magnetic field and trajactory evaluation as core of the code
libgtracr = Extension('gtracr.lib.libgtracr',
                      sources=[
                          "gtracr/lib/src/TrajectoryTracer.cpp",
                          "gtracr/lib/src/uTrajectoryTracer.cpp",
                          "gtracr/lib/src/igrf.cpp",
                          "gtracr/lib/src/pybind11_wrapper.cpp"
                      ],
                      language='c++',
                      include_dirs=['gtracr/lib/include'])
'''
Adopted from boost-histogram from scikit-HEP:
https://github.com/scikit-hep/boost-histogram/blob/develop/setup.py
'''
try:
    from numpy.distutils.ccompiler import CCompiler_compile
    import distutils.ccompiler

    distutils.ccompiler.CCompiler.compile = CCompiler_compile
except ImportError:
    print("Numpy not found, parallel compile not available")

# ext_modules = [
#     Extension("boost_histogram._core",
#               SRC_FILES,
#               include_dirs=INCLUDE_DIRS,
#               language="c++")
# ]


# As of Python 3.6, CCompiler has a `has_flag` method.
# cf http://bugs.python.org/issue26689
def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile

    with tempfile.NamedTemporaryFile("w", suffix=".cpp") as f:
        f.write("int main (int argc, char **argv) { return 0; }")
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except setuptools.distutils.errors.CompileError:
            return False
    return True


def cpp_flag(compiler):
    """Return the -std=c++14 compiler flag.
    """
    if has_flag(compiler, "-std=c++11"):
        return "-std=c++11"
    else:
        raise RuntimeError(
            "Unsupported compiler -- at least C++11 support is needed!")


class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""

    c_opts = {"msvc": ["/EHsc"], "unix": ["-g0"]}
    # c_opts = {"unix": ["-g0"]}

    if sys.platform == "darwin":
        c_opts["unix"] += ["-stdlib=libc++", "-mmacosx-version-min=10.9"]

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        # if ct == "unix":
        #     opts.append('-DVERSION_INFO="%s"' %
        #                 self.distribution.get_version())
        #     opts.append(cpp_flag(self.compiler))
        #     if has_flag(self.compiler, "-fvisibility=hidden"):
        #         opts.append("-fvisibility=hidden")
        # elif ct == "msvc":
        #     opts.append('/DVERSION_INFO=\\"%s\\"' %
        #                 self.distribution.get_version())
        for ext in self.extensions:
            ext.extra_compile_args = opts
        build_ext.build_extensions(self)


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

extras_require = {"examples": ["matplotlib", "mpld3"]}

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
    cmdclass={"build_ext": BuildExt},
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
