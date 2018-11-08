#!/usr/bin/env python

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(name='hugs',
      version='0.1',
      author='Johnny Greco',
      packages=['hugs'],
      url='https://github.com/johnnygreco/hugs',
      description='Hunt for Ultra-lsb Galaxies')
