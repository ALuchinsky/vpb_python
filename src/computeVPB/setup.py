from setuptools import setup
from Cython.Build import cythonize

setup(
    ext_modules=cythonize([
         'computeVPB.pyx'
         ],  # Python code file with primes() function
        annotate=True),                 # enables generation of the html annotation file
)