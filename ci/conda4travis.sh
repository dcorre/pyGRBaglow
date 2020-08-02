#!/bin/bash

# begin installing miniconda
if [[ "$TRAVIS_OS_NAME" != "windows" ]]; then
    echo "installing miniconda for posix";
    bash $HOME/download/miniconda.sh -b -u -p $MINICONDA_PATH;
elif  [[ "$TRAVIS_OS_NAME" == "windows" ]]; then
    echo "folder $MINICONDA_SUB_PATH does not exist"
    echo "installing miniconda for windows";
    if [[ ${PYTHON_VERSION} < 3 ]]; then
	    miniconda_version="miniconda";
    else
            miniconda_version="miniconda3";
    fi;
    echo "$miniconda_version";
    choco install $miniconda_version --params="'/JustMe /AddToPath:1 /D:$MINICONDA_PATH_WIN'";
fi;
# end installing miniconda

export PATH="$MINICONDA_PATH:$MINICONDA_SUB_PATH:$MINICONDA_LIB_BIN_PATH:$PATH";

# begin checking miniconda existance
echo "checking if folder $MINICONDA_SUB_PATH exists"
if [[ -d $MINICONDA_SUB_PATH ]]; then
    echo "folder $MINICONDA_SUB_PATH exists"
else
    echo "folder $MINICONDA_SUB_PATH does not exist"
fi;
# end checking miniconda existance

source $MINICONDA_PATH/etc/profile.d/conda.sh;

hash -r;
echo $TRAVIS_OS_NAME
echo $CONDA_PYTHON
python --version

# get CONDA base path
CONDA_PATH=$(conda info --base)

# configure miniconda
conda config --set always_yes yes --set changeps1 no;
conda config --add channels conda-forge;
conda update --quiet --yes conda;

# Useful for debugging any issues with conda
conda info -a

# create environment for tests (if needed)
if [ ! -f ${CONDA_PATH}/envs/tets-env/conda-meta/history ]; then
    conda create --name test-env python=${PYTHON_VERSION} pip setuptools
fi

# install conda dependencies (based on pip requirements file)
conda run --name test-env \
#conda install --name test-env --quiet --yes --file requirements_dev.txt --update-all
conda install --name test-env --yes --file requirements_dev.txt --update-all

# activate the environment
. ${CONDA_PATH}/etc/profile.d/conda.sh
conda activate test-env

# install other dependencies
pip install hjson
