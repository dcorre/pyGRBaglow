#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

import sys
from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize

if sys.platform.startswith("win"):
    # windows
    define_macros = [('USE_OPENMP', None)]
    libraries = ['kernel32']
    extra_compile_args = ['/EHsc', '/O2', '/Gy', '/openmp']
    extra_link_args = None
elif sys.platform.startswith("darwin"):
    # mac osx
    define_macros = [('USE_OPENMP', None)]
    libraries = None
    #extra_compile_args = ['-lpthread', '-stdlib=libc++',
    #                      '-mmacosx-version-min=10.7', '-Xpreprocessor',
    #                      '-fopenmp']
    extra_compile_args = ['-lpthread', '-fopenmp']
    extra_link_args = ["-lomp"]
else:
    # linux
    define_macros = [('USE_OPENMP', None)]
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
        'matplotlib']

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest', ]

extensions = [ 
        Extension("pyGRBaglow.synchrotron_model", ["pyGRBaglow/synchrotron_model.pyx"],
                  libraries=libraries,
                  extra_compile_args=extra_compile_args,
                  extra_link_args=extra_link_args)
        ]
setup(
    author="David Corre, Alain Klotz",
    author_email='david.corre.fr@gmail.com',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    description="GRB afterglow modelling with standard synchrotron model",
    entry_points={
        'console_scripts': [
            'pyGRBaglow=pyGRBaglow.cli:main',
        ],
    },
    ext_modules = cythonize(extensions),
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='pyGRBaglow',
    name='pyGRBaglow',
    packages=find_packages(),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/dcorre/pyGRBaglow',
    version='0.1.0',
    zip_safe=False,
)
