"""Test the reddening module."""

from extinctions import reddening


# Test values
RA = 340.83
DEC = -9.59


# Test the Reddening class

def test_init_float():
    """Initialize with two floats."""
    reddening.Reddening(RA, DEC)


def test_init_list():
    """Initialize with two list."""
    reddening.Reddening([RA], [DEC])


def test_init_mixt():
    """Initialize with one float and one list."""
    reddening.Reddening(RA, [DEC])


def test_astroquery():
    """Test the astroquery wrapper."""
    red = reddening.Reddening(RA, DEC)
    red.from_astroquery(dustmap='SFD98')
    red.from_astroquery(dustmap='SF11')


def test_snfactory():
    """Test the 'snfactory' query method."""
    red = reddening.Reddening(RA, DEC)
    red.from_snfactory()


def test_argonaut():
    """Test the query to the argonaut server."""
    red = reddening.Reddening(RA, DEC)
    red.from_argonaut()


def test_sncosmo():
    """Test the sncosmo wrapper."""
    red = reddening.Reddening(RA, DEC)
    red.from_sncosmo()


def test_sfd_map():
    """Test the local query on the SFD map."""
    red = reddening.Reddening(RA, DEC)
    red.query_local_map(dustmap='sfd')


def test_schlafly_map():
    """Test the local query on the Schlafly map."""
    red = reddening.Reddening(RA, DEC)
    red.query_local_map(dustmap='schlafly')
