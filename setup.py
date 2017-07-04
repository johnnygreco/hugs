#!/usr/bin/env python

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(name='hugs',
      version='v0.1',
      author='Johnny Greco',
      author_email='jgreco@astro.princeton.edu',
      packages=['hugs'],
      url='https://github.com/johnnygreco/hugs',
      description='Hunt for Ultra-lsb Galaxies')
