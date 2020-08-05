.. highlight:: shell

============
Installation
============

You can install it with pip:

.. code-block:: console

    pip install pyGRBaglow

It requires a compiler because of the usage of cython and openmp.

On MacOs, if it fails check if you have installed llvm and libomp. Otherwise intall them:

.. code-block:: console

    brew install llvm libomp

    export CPP=/usr/local/opt/llvm/bin/clang;


For Windows, only thing I know is that it builds successfully in Travis. So if installation using pip fails, look at the setup.py to have information about the compiler flags and build it yourself.


Installation from sources
-------------------------

The sources for pyGRBaglow can be downloaded from the `Github repo`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone git://github.com/dcorre/pyGRBaglow

Or download the `tarball`_:

.. code-block:: console

    $ curl -OL https://github.com/dcorre/pyGRBaglow/tarball/master

You can now install it with:

.. code-block:: console

    $ python setup.py [develop | install]

*develop* if you want to add modifications to the code. *install* otherwise.

.. _Github repo: https://github.com/dcorre/pyGRBaglow
.. _tarball: https://github.com/dcorre/pyGRBaglow/tarball/master


