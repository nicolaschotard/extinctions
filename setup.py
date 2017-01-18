#!/usr/bin/env python

"""Setup script."""

import os
import glob
import yaml
from setuptools import setup, find_packages


README = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/README.rst'
VERSION = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/version.yaml'

# Get version from version.py without importing package itself.
VERSION = yaml.load(open(VERSION))['version']

# Package name
NAME = 'extinctions'

# Packages (subdirectories in extinction/)
PACKAGES = find_packages()

# Scripts (in scripts/)
SCRIPTS = glob.glob("scripts/*.py")

PACKAGE_DATA = {NAME: ["data/maps.yaml"]}

setup(name=NAME,
      version=VERSION,
      description=("Extinction laws, maps and corrections"),
      license="MIT",
      classifiers=["Topic :: Scientific :: Astronomy",
                   "Intended Audience :: Science/Research"],
      url="https://github.com/nicolaschotard/extinctions",
      author="Nicolas Chotard",
      author_email="nchotard@in2p3.fr",
      packages=PACKAGES,
      scripts=SCRIPTS,
      package_data=PACKAGE_DATA,
      long_description=open(README).read(),
      setup_requires=['pytest-runner'],
      tests_require=['pytest']
     )
