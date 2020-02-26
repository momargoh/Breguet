import os
from setuptools import setup, find_packages
from setuptools.extension import Extension
import numpy

import Cython
from Cython.Distutils import build_ext
from Cython.Build import cythonize

cmdclass = {}
ext_modules = []
include_dirs = [numpy.get_include()]
compiler_directives = {
    'embedsignature': True,
    "language_level": 3,
    #"boundscheck": False,
    #"wraparound": False
}


def scanForExtension(directory, extension, files=[]):
    "Find all files with extension in directory and any subdirectories, modified from https://github.com/cython/cython/wiki/PackageHierarchy"
    length = len(extension)
    for f in os.listdir(directory):
        path = os.path.join(directory, f)
        if os.path.isfile(path) and path.endswith(extension):
            files.append(path[:-length])
        elif os.path.isdir(path):
            scanForExtension(path, extension, files)
    return files


pyx_exts = scanForExtension(".", ".pyx")
for ext in pyx_exts:
    ext_modules += cythonize(
        "{}.pyx".format(ext), compiler_directives=compiler_directives)
cmdclass.update({'build_ext': build_ext})

setup(
    name='breguet',
    packages=find_packages(),
    install_requires=['numpy', 'matplotlib', 'Cython'],
    extras_require={},
    python_requires='>=3',
    data_files=[],
    entry_points={},
    cmdclass=cmdclass,
    include_dirs=include_dirs,
    ext_modules=ext_modules,
    zip_safe=False,
)
