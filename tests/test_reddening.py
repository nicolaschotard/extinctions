"""Test the reddening module."""

import pytest
from Extinction import reddening


RA = 340.83
DEC = -9.59


@pytest.fixture()
def float_coordinates():
    return RA, DEC


@pytest.fixture()
def list_coordinates():    
    return [RA], [DEC]


@pytest.fixture()
def mixt_coordinates():    
    return RA, [DEC]


@pytest.fixture()
def get_reddening(float_coordinates):
    ra, dec = float_coordinates
    return reddening.Reddening(ra, dec)


class TestReddening:
    def test_init_float(self, float_coordinates):
        ra, dec = float_coordinates
        reddening.Reddening(ra, dec)

    def test_init_list(self, list_coordinates):
        ra, dec = list_coordinates
        reddening.Reddening(ra, dec)

    def test_init_mixt(self, mixt_coordinates):
        ra, dec = mixt_coordinates
        reddening.Reddening(ra, dec)

    def test_astroquery(self, get_reddening):
        get_reddening.from_astroquery(dustmap='SFD98')
        get_reddening.from_astroquery(dustmap='SF11')

    def test_snfactory(self, get_reddening):
        get_reddening.from_snfactory()

    def test_argonaut(self, get_reddening):
        get_reddening.from_argonaut()

    def test_sncosmo(self, get_reddening):
        get_reddening.from_sncosmo()

    def test_sfd_map(self, get_reddening):
        get_reddening.query_local_map(dustmap='sfd')

    def test_sncosmo(self, get_reddening):
        get_reddening.query_local_map(dustmap='schlafly')

    
