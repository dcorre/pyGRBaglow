#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `pyGRBaglow` package."""

import pytest
import numpy as np
import numpy.testing as npt

from click.testing import CliRunner
from pyGRBaglow import cli
from pyGRBaglow.reddening import reddening


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


def test_synchroton_1():
    """ Test a single simple run  """
    from pyGRBaglow.synchrotron_model import fireball_afterglow as grb
    import pyGRBaglow.constants as cc

    # define GRB parameters
    redshift = 3.92
    n0 = 258.08
    eps_b = 0.0272
    eps_e = 0.547
    E_iso = 4.41e53
    eta = 0.77
    p = 2.68
    Y = 0.
    ism_type = 0.
    # in angstroms
    wavelength = np.logspace(-7, 12, 10000)
    # in days
    time = np.array([30/86400, 10/1440, 1/24, 1])
    frequencies = 3e8 / (wavelength * 1e-10)

    # Load object
    afterglow = grb(n0=n0, eps_b=eps_b, eps_e=eps_e, E_iso=E_iso,
                    eta=eta, p=p, Y=Y, z=redshift, ism_type=ism_type,
                    disp=0, num_threads=4)
    # Compute light curve for each time
    afterglow_lc = afterglow.light_curve(time, frequencies)

    expected_lc = np.genfromtxt('tests/data/lc_synchrotron.dat')
    npt.assert_allclose(expected_lc, afterglow_lc.T, rtol=1e-6, atol=0)


def test_synchrotron_2():
    """
    Test a single simple run with windy environment and inverse compton
    """
    from pyGRBaglow.synchrotron_model import fireball_afterglow as grb
    import pyGRBaglow.constants as cc

    # define GRB parameters
    redshift = 1
    n0 = 1e-2
    eps_b = 1e-5
    eps_e = 1e-3
    E_iso = 4.41e50
    eta = 0.8
    p = 2.2
    Y = 0.5
    ism_type = 2
    # in angstroms
    wavelength = np.logspace(-7, 12, 10000)
    # in days
    time = np.array([30/86400, 1/24, 1, 10])
    frequencies = 3e8 / (wavelength * 1e-10)

    # Load object
    afterglow = grb(n0=n0, eps_b=eps_b, eps_e=eps_e, E_iso=E_iso,
                    eta=eta, p=p, Y=Y, z=redshift, ism_type=ism_type,
                    disp=0, num_threads=4)
    # Compute light curve for each time
    afterglow_lc = afterglow.light_curve(time, frequencies)

    expected_lc = np.genfromtxt('tests/data/lc_synchrotron_windy.dat')
    npt.assert_allclose(expected_lc, afterglow_lc.T, rtol=1e-6, atol=0)


@pytest.mark.parametrize("template", ['SPL', 'BPL'])
def test_template(template):
    """ Test a single simple run with SPL and BPL """
    from pyGRBaglow.template_models import Templates

    # define GRB parameters
    # in angstroms
    wavelength = np.logspace(-7, 12, 10000)
    # in days
    time = np.array([30/86400, 1/24, 1, 10])

    if template == 'SPL':
        params = [-1, 1.5]
    elif template == 'BPL':
        params = [-1, -0.5, 1.5, 1]

    # Load object
    temp = Templates()
    lc = temp.light_curve(wavelength, time, params, model=template)

    expected_lc = np.genfromtxt('tests/data/test_lc_%s.dat' % template)
    npt.assert_allclose(expected_lc, lc.T, rtol=1e-6, atol=0)


@pytest.mark.parametrize("extLaw", ['mw', 'lmc', 'smc'])
def test_dust_extinction_Pei(extLaw):

    wavelength = np.logspace(-7, 12, 10000)
    redshift = 1.5
    Av = 0.5
    red_model = reddening(wavelength, redshift, Av)
    trans_dust_host = red_model.Pei92(ext_law=extLaw, Xcut=True)[1]
    expected_trans = np.genfromtxt('tests/data/test_reddening_Pei_%s.dat'
                                   % extLaw)
    npt.assert_allclose(expected_trans, trans_dust_host)


@pytest.mark.parametrize("extLaw", ['mw', 'lmc', 'smc',
                                    'linear', 'calzetti', 'grb1', 'grb2'])
def test_dust_extinction_Li(extLaw):

    wavelength = np.logspace(-7, 12, 10000)
    redshift = 1.5
    Av = 0.5
    trans_dust_host = reddening(wavelength, redshift, Av).Li07(extLaw,
                                                               Xcut=True)[1]
    expected_trans = np.genfromtxt('tests/data/test_reddening_Li_%s.dat'
                                   % extLaw)
    npt.assert_almost_equal(expected_trans, trans_dust_host, decimal=3)


@pytest.mark.parametrize("NHx", [0.2, 1])
def test_gas_absorption(NHx):

    wavelength = np.logspace(-7, 12, 10000)
    redshift = 1.5
    Av = 0.5
    npt.assert_allclose(4 * Av * 1.79e21,
                        reddening(wavelength, 1.5, Av=0.5).NHI_host())

    trans = reddening(wavelength, 1.5, Av=0.5).gas_absorption(NHx=NHx)
    if NHx == 0.2:
        label = '02'
    elif NHx == 1:
        label = '1'
    expected_trans = np.genfromtxt('tests/data/test_gas_abs_%s.dat' % label)
    npt.assert_allclose(expected_trans, trans)


@pytest.mark.parametrize("z", [0.5, 2, 10])
@pytest.mark.parametrize("igm_model", ['meiksin', 'madau'])
def test_igm(igm_model, z):
    from pyGRBaglow.igm import meiksin, madau
    wavelength = np.logspace(-7, 12, 10000)
    if igm_model == 'meiksin':
        trans = meiksin(wavelength/10, z, Xcut=True)
    elif igm_model == 'madau':
        trans = madau(wavelength, z, lylim=True, metals=True, Xcut=True)
    if z == 0.5:
        label = '05'
    elif z == 2:
        label = '2'
    elif z == 10:
        label = '10'
    expected_trans = np.genfromtxt('tests/data/test_igm_%s_%s.dat' %
                                   (igm_model, label))
    npt.assert_allclose(expected_trans, trans)


@pytest.mark.parametrize("z", [1.5, 5])
@pytest.mark.parametrize("NHI", [1e15, 1e20])
def test_dla(z, NHI):
    from pyGRBaglow.igm import dla
    wavelength = np.logspace(-7, 12, 10000)
    trans = dla(wavelength, z, NHI)
    if z == 1.5:
        label = '15'
    elif z == 5:
        label = '5'
    if NHI == 1e15:
        label2 = '1e15'
    elif NHI == 1e20:
        label2 = '1e20'
    expected_trans = np.genfromtxt('tests/data/test_dla_%s_%s.dat' %
                                   (label, label2))
    npt.assert_allclose(expected_trans, trans)
