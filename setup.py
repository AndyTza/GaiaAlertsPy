#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os

from setuptools import setup

setup(
    name='GaiaAlertsPy',
    version='0.1.1',
    author='Andy Tzanidakis',
    author_email='atzanida@uw.edu',
    description='Fetch Gaia Photometric Science Alerts',
    long_description='Fetch Gaia Photometric Science Alerts',
    url='https://github.com/AndyTza/GaiaAlertsPy',
    license='BSD-3-Clause',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    packages=['GaiaAlertsPy']
)
