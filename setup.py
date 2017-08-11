#!/usr/bin/env python
from setuptools import setup, find_packages
__author__ = 'adamkoziol'

setup(
    name="antismash-parse",
    version="0.1",
    packages=find_packages(),
    include_package_data=True,
    license='MIT',
    author='Adam Koziol',
    author_email='adam.koziol@inspection.gc.ca',
    description='antiSMASH parsing software',
    url='https://github.com/adamkoziol/antismash-parse.git',
    long_description=open('README.md').read(),
    install_requires=['biopython >= 1.65',
                      'OLCTools'],
)
