#!/usr/bin/env python

from setuptools import setup, find_packages


def readfile(filename):
    with open(filename, 'r+') as f:
        return f.read()


setup(name='gammagwfollowup',
      version=0.1,
      description="Gravitational waves follow-up scheduling and simulation for IACTs",
      install_requires=[
          'gammapy==0.9',
          'astropy==3.0.5',
          'scipy',
          'healpy==1.12.9',
          'ipython',
          'matplotlib==2.2.2',
          'MOCpy==0.6.0',
          'numpy',
          'pandas',
          'pytz==2021.1',
          'ephem==3.7.6.0'
      ],
      packages=find_packages(),
      # tests_require=['pytest'],
      author='Monica Seglar-Arroyo ',
      author_email='mosear2@gmail.com ',
      url='https://drf-gitlab.cea.fr/multimessenger-IRFU/cta/gw-follow-up-simulations',
      long_description=readfile('README.md'),
      include_package_data=True,
      package_data={'': ['dataset/*.all']},
      # license='',
      classifiers=[],
      scripts=[],
      )
