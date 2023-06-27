#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os

from setuptools import setup

VERSION_TEMPLATE = """
try:
    from setuptools_scm import get_version
    __version__ = get_version(root='..', relative_to=__file__)
except Exception:
    __version__ = '1.0.0'
""".lstrip()

setup(
    use_scm_version={'write_to': os.path.join('GaiaAlertsPy', 'history.py'),
                     'write_to_template': VERSION_TEMPLATE},
                         name='GaiaAlertsPy',
                        version='1.0.0',
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
                        packages=['GaiaAlertsPy'])