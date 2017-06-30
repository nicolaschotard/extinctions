"""Test the scripts."""

import os
import matplotlib as mpl
if os.environ.get('DISPLAY', '') == '':
    print 'no display found. Using non-interactive Agg backend'
    mpl.use('Agg')
from extinctions import main


def test_get_maps():
    """Test get_maps.py."""
    main.get_maps(["--list"])
    main.get_maps(["--update", "--select", "sfd"])
    main.get_maps(["--update", "--exclude", "green,planck,schlafly,sfd_npg,sfd_spg"])


def test_extinction_plots():
    """Test extinction_plots.y."""
    main.extinction_plots(["--hide"])
