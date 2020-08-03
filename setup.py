#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize

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
                  libraries=["m"],
                  extra_compile_args=["-ffast-math", "-fopenmp"],
                  extra_link_args=["-ffast-math", "-fopenmp"])
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
