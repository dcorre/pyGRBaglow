#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `pyGRBaglow` package."""

import pytest
import numpy.testing as npt

from click.testing import CliRunner

from pyGRBaglow import pyGRBaglow
from pyGRBaglow import cli


@pytest.fixture
def response():
    """Sample pytest fixture.

    See more at: http://doc.pytest.org/en/latest/fixture.html
    """
    # import requests
    # return requests.get('https://github.com/audreyr/cookiecutter-pypackage')


def test_content(response):
    """Sample pytest test function with the pytest fixture as an argument."""
    # from bs4 import BeautifulSoup
    # assert 'GitHub' in BeautifulSoup(response.content).title.string


def test_command_line_interface():
    """Test the CLI."""
    runner = CliRunner()
    result = runner.invoke(cli.main)
    assert result.exit_code == 0
    assert 'pyGRBaglow.cli.main' in result.output
    help_result = runner.invoke(cli.main, ['--help'])
    assert help_result.exit_code == 0
    assert '--help  Show this message and exit.' in help_result.output

def test_simple_run()
    """ Test a single simple run  """
    from pyGRBaglow.synchrotron_model import fireball_afterglow as grb
    import pyGRBaglow.constants as cc

    # define GRB parameters
    redshift = 3.92
    n0 = 258.08
    eps_b = 0.0272
    eps_e = 0.547
    E_iso = 4.41e53
    eta=0.77
    p=2.68  #>2
    Y=0
    ism_type=0

    wavelength=np.logspace(-7,12,1e4) #in angstroms
    time = [30/86400,10/1440,1/24,1] # in days
    frequencies = 3e8 / (wavelength*1e-10)

    #Load object
    afterglow=grb(n0=n0, eps_b=eps_b, eps_e=eps_e, E_iso=E_iso, eta=eta, p=p,
                  Y=Y, z=redshift, ism_type=ism_type, disp=0)
    #Compute light curve for each time
    afterglow_lc=afterglow.light_curve(time, frequencies)

    expected_lc = np.genfromtxt('lc_test')

    npt.assert_equal(expected_lc, afterglow_lc.T)

