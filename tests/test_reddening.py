"""Test the reddening module."""

import pytest
from Extinction import reddening


RA = 340.83
DEC = -9.59


@pytest.fixture()
def float_coordinates():
    """Return RA DEC as two floats."""
    return RA, DEC


@pytest.fixture()
def list_coordinates():
    """Return RA DEC as two lists."""
    return [RA], [DEC]


@pytest.fixture()
def mixt_coordinates():
    """Return RA as a float, DEC as a list."""
    return RA, [DEC]


@pytest.fixture()
def get_reddening(float_coordinates):
    """Get the Reddening object."""
    ra, dec = float_coordinates
    return reddening.Reddening(ra, dec)


class TestReddening(object):

    """Test the Reddening class."""

    def __init__(self):
        print "Testing reddening.Reddening"
        
    def test_init_float(self, float_coordinates):
        """Initialize with two floats."""
        ra, dec = float_coordinates
        reddening.Reddening(ra, dec)

    def test_init_list(self, list_coordinates):
        """Initialize with two list."""
        ra, dec = list_coordinates
        reddening.Reddening(ra, dec)

    def test_init_mixt(self, mixt_coordinates):
        """Initialize with one float and one list."""
        ra, dec = mixt_coordinates
        reddening.Reddening(ra, dec)

    def test_astroquery(self, get_reddening):
        """Test the astroquery wrapper."""
        get_reddening.from_astroquery(dustmap='SFD98')
        get_reddening.from_astroquery(dustmap='SF11')

    def test_snfactory(self, get_reddening):
        """Test the 'snfactory' query method."""
        get_reddening.from_snfactory()

    def test_argonaut(self, get_reddening):
        """Test the query to the argonaut server."""
        get_reddening.from_argonaut()

    def test_sncosmo(self, get_reddening):
        """Test the sncosmo wrapper."""
        get_reddening.from_sncosmo()

    def test_sfd_map(self, get_reddening):
        """Test the local query on the SFD map."""
        get_reddening.query_local_map(dustmap='sfd')

    def test_schlafly_map(self, get_reddening):
        """Test the local query on the Schlafly map."""
        get_reddening.query_local_map(dustmap='schlafly')
