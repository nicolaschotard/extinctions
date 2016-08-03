import os
import pyfits
import healpy
import numpy
import astropy.units as units
from astropy.coordinates import SkyCoord

from ToolBox.Astro import Fetchers
from astroquery.irsa_dust import IrsaDust

def test():
    print "test"

def load_map(map=0):
    """
    map is eaither 0 (sfd) or 1 (schlafly)
    """
    if map == 0:
        map = os.getenv('HOME') + '/.extinction/maps/lambda_sfd_ebv.fits'
        hmap = healpy.read_map(map)
    elif map == 1:
        map = os.getenv('HOME') + '/.extinction/maps/ps1-ebv-4.5kpc.fits'
        map = pyfits.getdata(map)
        hmap = map['ebv']
    return hmap

def plot_map(map=0):
    """
    map is eaither 0 (sfd) or 1 (schlafly)
    """
    if map == 0:
        title = 'SFD E(B-V) map'
    elif map == 1:
        title = 'Schlafly E(B-V) map'

    map = load_map(map)
    healpy.mollview(map, title=title, unit='mag',
                    norm='hist', min=0, max=0.5, xsize=2000)
    healpy.graticule()
    #healpy.gnomview(map, rot=[0,0.3], title='GnomView', unit='mK', format='%.2g')

def get_value(ra, dec, map=0):

    nest = 1 if map else 0
    
    # Parse input
    coordinates = SkyCoord(ra=ra, dec=dec,
                           unit=units.degree) 

    # Convert to galactic coordinates.
    l = coordinates.galactic.l.degree
    b = coordinates.galactic.b.degree
    map = load_map(map)
    ebv = healpy.get_interp_val(map,  (90.-b)*numpy.pi/180., l*numpy.pi/180.,
                                nest=nest)
    n = healpy.get_interp_weights(512, (90.-b)*numpy.pi/180., l*numpy.pi/180.,
                                  nest=nest)
    print n, 'toto'
    print map[5000]
    print numpy.sum([map[i]*j for i, j in zip(n[0], n[1])])
    ebv2 = numpy.mean([map[i] for i in n[0]])
    t = IrsaDust.get_extinction_table('%.4f %.4f' % (ra, dec))
    f = Fetchers.sfd_ebmv(ra, dec)
    print ebv, ebv2, t[9]['A_SFD'] / t[9]['A_over_E_B_V_SFD'], f
