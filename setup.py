import os
import glob
from distutils.core import setup

scripts=['meds-extract-range']

scripts=[os.path.join('./scripts', s) for s in scripts]

setup(name="meds", 
      version="0.1.0",
      description="Python and C libraries for reading MEDS files",
      license = "GPL",
      author="Erin Scott Sheldon",
      author_email="erin.sheldon@gmail.com",
      packages=['meds'],
      scripts=scripts)

