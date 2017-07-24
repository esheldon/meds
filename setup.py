import os
import glob
from distutils.core import setup, Extension
import numpy

scripts=[
    'meds-extract-range',
    'meds-extract-catalog',
    'meds-view',
]
scripts=[os.path.join('./scripts', s) for s in scripts]

sources=["meds/_uberseg.c"]
include_dirs=[numpy.get_include()]
ext=Extension("meds._uberseg", sources, include_dirs=include_dirs)

setup(name="meds", 
      version="0.9.3",
      description="Python and C libraries for reading MEDS files",
      license = "GNU GPLv3",
      author="Erin Scott Sheldon",
      author_email="erin.sheldon@gmail.com",
      packages=['meds'],
      ext_modules=[ext],
      include_dirs=include_dirs,
      scripts=scripts)

