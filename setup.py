from __future__ import print_function

import os
import sys
import subprocess
import shutil

try:
    import setuptools
except ImportError:
    sys.exit("setuptools package not found. "
             "Please use 'pip install setuptools' first")

from setuptools import setup

# Make sure we're running from the setup.py directory.
script_dir = os.path.dirname(os.path.realpath(__file__))
if script_dir != os.getcwd():
    os.chdir(script_dir)

from hapdup.__version__ import __version__


setup(name="hapdup",
      version=__version__,
      description="Pipeline to convert a haploid assembly into diploid",
      url="https://github.com/fenderglass/HapDup",
      author="Mikhail Kolmogorov",
      author_email = "fenderglass@gmail.com",
      license="BSD-3-Clause",
      packages=["hapdup"],
      entry_points={"console_scripts": ["hapdup = hapdup.main:main"]}
      )

