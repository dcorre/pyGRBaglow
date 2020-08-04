# pyGRBaglow: GRB afterglow modelling with standard synchrotron model 
  
* Free software: MIT license
* Documentation: https://pyGRBaglow.readthedocs.io.


Release status
--------------

[![PyPI version](https://badge.fury.io/py/pyGRBaglow.svg)](http://badge.fury.io/py/pyGRBaglow)
![Supported Python versions](https://img.shields.io/pypi/pyversions/pyGRBaglow.svg)


Development status
--------------------

[![Build Status](https://travis-ci.com/dcorre/pyGRBaglow.svg?branch=master)](https://travis-ci.com/dcorre/pyGRBaglow)
[![codecov](https://codecov.io/gh/dcorre/pyGRBaglow/branch/master/graphs/badge.svg)](https://codecov.io/gh/dcorre/pyGRBaglow/branch/master)
[![Documentation Status](https://readthedocs.org/projects/pygrbaglow/badge/?version=latest)](https://pyetc.readthedocs.io/en/latest/?badge=latest)

Features
--------
* Compute Gamma-Ray Burst light curves using:
    * standard synchrotron model of [Granot&Sari+02](https://ui.adsabs.harvard.edu/abs/2002ApJ...568..820G/abstract).
    * empirical templates (Single power-law, Broken power-law)

* These simulated light curves are extinguished along the line of sight due to the type and amount of dust in the host galaxy, and due to photons absorbed in the Intergalactic medium. 
    * Different extinction laws from [Pei+92](http://adsabs.harvard.edu/abs/1992ApJ...395..130P) to compute amount of extinction along one line of sight occuring for a given wavelength at a given redshift.
    * Intergalactic medium (IGM) mean transmission along one line of sight for a given redshift for a given redshift. 2 models available: [Madau+95](http://adsabs.harvard.edu/abs/1995ApJ...441...18M) and [Meiksin+06](http://adsabs.harvard.edu/abs/2006MNRAS.365..807M).
    * Photoelectric absorption model based on cross sections of [McCammon+93](http://adsabs.harvard.edu/abs/1983ApJ...270..119M).

* Use of cython to speed up the execution time, as well as parallelisation using openmp.

* Build successfully on Linux, macOS and Windows.


Installation
------------

pip install pyGRBaglow


Credits
-------

This package was created with [Cookiecutter](https://github.com/audreyr/cookiecutter) and the [audreyr/cookiecutter-pypackage](https://github.com/audreyr/cookiecutter-pypackage) project template.
 
