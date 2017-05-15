"""Get the reddening E(B-V) for a given set of coordinates."""

import os
import numpy as np
from astropy.io import fits
from astropy import units
from astropy.coordinates import SkyCoord
import healpy
import yaml
from pkg_resources import resource_filename

from extinctions.extern import snfactory
from extinctions.extern import argonaut
from extinctions.extern import sncosmo
from extinctions.extern import others


class Reddening(object):

    """Query reddening from different sources."""

    def __init__(self, ra, dec, map_dir=None, loadmaps=False):
        """Input are the rad/dec coordinates in degree and a map directory."""
        self.ra = ra if isinstance(ra, list) else [ra]
        self.dec = dec if isinstance(dec, list) else [dec]
        assert len(self.ra) == len(self.dec)
        self.coordinates = SkyCoord(ra=units.Quantity(ra, 'deg'),
                                    dec=units.Quantity(dec, 'deg'),
                                    unit=units.degree)

        # Convert to galactic coordinates.
        self.theta = (90. - self.coordinates.galactic.b.degree) * np.pi / 180.
        self.phi = self.coordinates.galactic.l.degree * np.pi / 180.

        # Map directory
        map_dir = os.environ.get("MAPSDIR") if map_dir is None else map_dir
        map_dir = os.getenv('HOME') + '/.extinction/maps/' if map_dir is None else map_dir
        self.map_dir = map_dir

        # Load the maps
        self.maps = self.loaded_maps = {}
        self._load_all_maps(loadall=loadmaps)

    def _load_all_maps(self, loadall=False):
        """Load the local maps."""
        pmap = self.map_dir + '/maps.yaml' if os.path.exists(self.map_dir + '/maps.yaml') else None
        pmap = resource_filename('extinctions', 'data/maps.yaml') if pmap is None else pmap
        if pmap is None:
            raise IOError("No maps.yaml found anywehre")
        print "INFO: Loading the maps from", pmap
        self.maps = yaml.load(open(pmap))
        if loadall:
            for m in self.maps.keys():
                self._load_maps(cmap=m)

    def _load_maps(self, cmap=None):
        lmap = self.map_dir + '/' + os.path.basename(self.maps[cmap]['url'])
        if not os.path.exists(lmap):
            print " - WARNING: You must download the map %s (%s) in order " % \
                (os.path.basename(lmap), cmap) + "to use it. Use get_maps to do so."
            self.maps.pop(cmap)
        elif cmap in ['sfd', 'planck']:
            field = 2 if cmap == 'planck' else 0
            self.loaded_maps[cmap] = self.maps[cmap]
            self.loaded_maps[cmap]['map'] = healpy.read_map(lmap, verbose=False,
                                                            field=field, memmap=True)
            print ' - ', cmap, "is loaded"
        elif cmap in ['schlafly', 'green']:
            self.loaded_maps[cmap] = self.maps[cmap]
            self.loaded_maps[cmap]['map'] = fits.getdata(lmap)['ebv']
            print ' - ', cmap, "is loaded"

    def from_astroquery(self, dustmap='SFD98'):
        """Query IRAS using the astropy/astroquery tools (SFD98 or SF11 maps)."""
        if len(self.ra) >= 2:
            print "WARNING: This online query is SLOW for several set of coordinates"
        tables = [others.astroquery.get_extinction_table('%.4f %.4f' % (ra, dec))
                  for ra, dec in zip(self.ra, self.dec)]
        if dustmap == 'SFD98':
            return [np.median(t['A_SFD'] / t['A_over_E_B_V_SFD']) for t in tables]
        else:
            return [np.median(t['A_SandF'] / t['A_over_E_B_V_SandF']) for t in tables]

    def from_snfactory(self):
        """Query on IRSA or NED using code from the SNfactory ToolBox."""
        if len(self.ra) >= 2:
            print "WARNING: This online query is SLOW for several set of coordinates"
        return [snfactory.sfd_ebmv(ra, dec) for ra, dec in zip(self.ra, self.dec)]

    def from_sncosmo(self):
        """Using sncosmo query utilisies on local SFD98 north/south maps."""
        return sncosmo.get_ebv_from_map([self.ra, self.dec], mapdir=self.map_dir)

    def from_argonaut(self):
        """Using the distant argonaut query utility."""
        return argonaut.query(self.ra, self.dec, coordsys='equ', mode='sfd')['EBV_SFD']

    def query_local_map(self, dustmap='sfd'):
        """Query one of the local map."""
        nest = True if dustmap in ['green'] else False
        if dustmap not in self.loaded_maps:
            self._load_maps(dustmap)
        return healpy.get_interp_val(self.loaded_maps[dustmap]['map'],
                                     self.theta, self.phi, nest=nest)


def load_map(lmap=0):
    """
    Load a map.

    map is either
    0: sfd
    1: schlafly
    2: std_s
    3: sfd_n
    4: planck
    """
    if lmap == 0:
        lmap = os.getenv('HOME') + '/.extinction/maps/lambda_sfd_ebv.fits'
        hmap = healpy.read_map(lmap)
    elif lmap == 1:
        lmap = os.getenv('HOME') + '/.extinction/maps/ps1-ebv-4.5kpc.fits'
        lmap = fits.getdata(lmap)
        hmap = lmap['ebv']
    elif lmap == 2:
        lmap = os.getenv('HOME') + '/.extinction/maps/SFD_dust_4096_sgp.fits'
        hmap = fits.getdata(lmap)
    elif lmap == 3:
        lmap = os.getenv('HOME') + '/.extinction/maps/SFD_dust_4096_ngp.fits'
        hmap = fits.getdata(lmap)
    elif lmap == 4:
        lmap = os.getenv('HOME') + '/.extinction/maps/HFI_CompMap_ThermalDustModel_2048_R1.20.fits'
        hmap = healpy.read_map(lmap)
    elif lmap == 5:
        lmap = os.getenv('HOME') + '/.extinction/maps/lambda_green_dust_map_2d.fits'
        lmap = fits.getdata(lmap)
        hmap = lmap['ebv']
    return hmap


def plot_map(smap=0):
    """
    Plot a map.

    smap is either
    0: sfd
    1: schlafly
    2: std_s
    3: sfd_n
    """
    if smap == 0:
        title = 'SFD E(B-V) map'
    elif smap == 1:
        title = 'Schlafly E(B-V) map'

    smap = load_map(smap)
    healpy.mollview(smap, title=title, unit='mag',
                    norm='hist', min=0, max=0.5, xsize=2000)
    healpy.graticule()


def test_ebm(ra, dec, smap=0, nest=False):
    """Make some tests."""
    # Parse input
    coordinates = SkyCoord(ra=ra, dec=dec, unit=units.degree)

    # Convert to galactic coordinates.
    l = coordinates.galactic.l.degree
    b = coordinates.galactic.b.degree
    theta = (90. - b) * np.pi / 180.
    phi = l * np.pi / 180.
    print "l, b = %.3f, %.3f" % (l, b)
    print "theta, phi = %.3f, %.3f" % (theta, phi)
    m = load_map(smap)

    # from this code
    if smap == 5:
        nest = True
    ebv = healpy.get_interp_val(m, theta, phi, nest=nest)

    # from astroquery
    t = others.astroquery.get_extinction_table('%.4f %.4f' % (ra, dec))
    if smap in [0, 2, 3]:
        t = t[9]['A_SFD'] / t[9]['A_over_E_B_V_SFD']
    else:
        t = t[9]['A_SandF'] / t[9]['A_over_E_B_V_SandF']
        print t

    # from SNf code (ned)
    f = snfactory.sfd_ebmv(ra, dec)

    # from sncosmo
    sn = sncosmo.get_ebv_from_map([ra, dec], mapdir='/home/chotard/.extinction/maps/')

    # from other query
    ebv_sfd = argonaut.query(ra, dec, coordsys='equ', mode='sfd')['EBV_SFD'][0]

    print "\nAll results:"
    print " - Healpy (lambda/nasa map): %.5f" % ebv
    print " - Astropy/IrsaDust: %.5f" % t
    print " - SNf code (irsa or ned): %.5f" % f, f
    print " - sncosmo (local N/S maps): %.5f" % sn
    print " - argonaut.skypams: %.5f" % ebv_sfd
