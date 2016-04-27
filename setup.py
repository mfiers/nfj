#!/usr/bin/env python

from setuptools import setup, find_packages

#one line description
with open('DESCRIPTION') as F:
    description = F.read().strip()

#version number
with open('VERSION') as F:
    version = F.read().strip()

entry_points = {
    'console_scripts': [
        'nfj = nfj.cli:dispatch'
        ]}

setup(name='nfj',
      version=version,
      description=description,
      author='mf',
      author_email='mf',
      entry_points = entry_points,
      include_package_data=True,
      url='https://encrypted.google.com/#q=nfj&safe=off',
      packages=find_packages(),
      install_requires=[
                'Leip',
                ],
      classifiers = [
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        ]
     )
