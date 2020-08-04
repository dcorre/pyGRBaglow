#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

import sys
from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize

# Lines adapted from https://github.com/sdpython/td3a_cpp
if sys.platform.startswith("win"):
    # windows
    libraries = ['kernel32']
    extra_compile_args = ['/EHsc', '/O2', '/Gy', '/openmp']
    extra_link_args = None
elif sys.platform.startswith("darwin"):
    # mac osx
    libraries = None
    extra_compile_args = ['-lpthread', '-stdlib=libc++',
                          '-mmacosx-version-min=10.7', '-Xpreprocessor',
                          '-fopenmp']
    extra_link_args = ["-lomp"]
else:
    # linux
    libraries = ["m"]
    extra_compile_args = ['-lpthread', '-fopenmp', "-ffast-math"]
    # option '-mavx2' forces the compiler to use
    # AVX instructions the processor might not have
    extra_link_args = ['-lgomp', "-ffast-math", "-fopenmp"]

with open('README.md') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
        'Click>=6.0',
        'numpy',
        'astropy',
        'scipy',
        'jupyter',
        'matplotlib',
        'cython']

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest', ]

extensions = [
        Extension("pyGRBaglow.synchrotron_model",
                  ["pyGRBaglow/synchrotron_model.pyx"],
                  libraries=libraries,
                  extra_compile_args=extra_compile_args,
                  extra_link_args=extra_link_args)
        ]
setup(
    author="David Corre",
    author_email='david.corre.fr@gmail.com',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="GRB afterglow modelling with standard synchrotron"
                "model of Granot & Sari 2002",
    entry_points={
        'console_scripts': [
            'pyGRBaglow=pyGRBaglow.cli:main',
        ],
    },
    ext_modules=cythonize(extensions),
    install_requires=requirements,
    license="MIT license",
    long_description=readme,
    long_description_content_type='text/markdown',
    include_package_data=False,
    keywords=['GRB', 'afterglow', 'fireball', 'synchrotron', 'light curve'],
    name='pyGRBaglow',
    packages=find_packages(),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/dcorre/pyGRBaglow',
    download_url='https://github.com/dcorre/pyGRBaglow/archive/v0.1.3.tar.gz',
    version='0.1.3',
    zip_safe=False,
)
