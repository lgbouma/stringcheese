# -*- coding: utf-8 -*-

"""
setup.py - Luke Bouma (luke@astro.princeton.edu) - Apr 2019

Stolen from the astrobase setup.py
"""
__version__ = '0.0.0'

import sys
from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

INSTALL_REQUIRES = [
    'numpy>=1.4.0',
    'scipy',
    'astropy>=1.3',
    'matplotlib',
]


###############
## RUN SETUP ##
###############

# run setup.
setup(
    name='stringcheese',
    version=__version__,
    description=('do strings of a feather spin together?'),
    classifiers=[
        'Development Status :: 1 - Alpha',
        'License :: OSI Approved :: MIT License',
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
    ],
    keywords='astronomy',
    url='https://github.com/lgbouma/stringcheese',
    author='Luke Bouma',
    author_email='luke@astro.princeton.edu',
    license='MIT',
    packages=[
        'stringcheese',
        'stringcheese.wrangling'
    ],
    install_requires=INSTALL_REQUIRES,
    include_package_data=False,
    zip_safe=False
)
