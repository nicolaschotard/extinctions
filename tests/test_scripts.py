"""Test the scripts."""

from extinctions import main


def test_get_maps():
    """Test get_maps.py."""
    main.get_maps(["--list"])
    main.get_maps(["--update", "--select", "sfd"])
    main.get_maps(["--update", "--exclude", "green,planck,schlafly,sfd_npg,sfd_spg"])


def test_extinction_plots():
    """Test extinction_plots.y."""
    main.extinction_plots(["--hide"])
