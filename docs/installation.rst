.. highlight:: shell

============
Installation
============


Dependencies
------------

You must have already installed the following packages on your computer:

- numpy

Create python environment (optional but advised)
------------------------------------------------

This will create locally a python environment in which you can install and specify the required libraries.
The advantage is that it will avoid you to change the libraries version you are using for other projects.
You need to install conda first:

- Download Miniconda for python 3 here: https://conda.io/miniconda.html (recommended)

- If you want to download a complete python environment (~1.5GB) you can download Anaconda (instead of Miniconda) for Python 3 here: https://www.continuum.io/downloads

Then to create the python environment, open a terminal and type:

.. code-block:: console

    conda create --name colibri python=3  numpy


We named this environment *colibri* as this the telescope for which we use this Image Simulator, but you can use an other environment name.

Once it is installed, type in a terminal (if you let *colibri* as the environment name, otherwise write the one you chose):

.. code-block:: console

    source activate colibri


It will activate the environment. You can see that (*colibri*) is added in front of your ID in your terminal. If you type "conda list" you can check which libraries and version are installed. When you want to exit this environment type "source deactivate".



Installation from sources
-------------------------

The sources for pyGRBaglow can be downloaded from the `Github repo`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone git://github.com/dcorre/pyGRBaglow

Or download the `tarball`_:

.. code-block:: console

    $ curl  -OL https://github.com/dcorre/pyGRBaglow/tarball/master


Once you have a copy of the source, if you created a python environment, do not forget to activate it with:

.. code-block:: console

    source activate colibri

Before installing it, remember that you need numpy to be installed.

You can now install it with:

.. code-block:: console

    $ python setup.py develop


.. _Github repo: https://github.com/dcorre/pyGRBaglow
.. _tarball: https://github.com/dcorre/pyGRBaglow/tarball/master


