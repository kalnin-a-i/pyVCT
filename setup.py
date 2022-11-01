from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import os
import shutil


extension = Extension(
    name="pyVCT",
    sources=["libcpmfem.pyx"],
    libraries=["cpmfem"],
    library_dirs=["VCT"],
    include_dirs=["VCT"]
)
setup(
    name="pyVCT",
    ext_modules=cythonize([extension])
)

os.remove('libcpmfem.c')
shutil.rmtree('build')