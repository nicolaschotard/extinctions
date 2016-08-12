"""
Set of utilities function
"""

import os
import healpy
import numpy as np

from astropy.io import fits
import astropy.units as units
from astropy.coordinates import SkyCoord

from Extinction.extern import snfactory, argonaut
from astroquery.irsa_dust import IrsaDust as astroquery
from sncosmo import dustmap as sncosmo

def load_map(lmap=0):
    """
    map is either 
    0: sfd
    1: schlafly
    2: std_s
    3: sfd_n
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
    return hmap

def plot_map(map=0):
    """
    Plot a map
    
    map is either
    0: sfd
    1: schlafly
    2: std_s
    3: sfd_n
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

def test_ebm(ra, dec, map=0):
    """
    Make some tests
    """
    # Parse input
    coordinates = SkyCoord(ra=ra, dec=dec,
                           unit=units.degree) 

    # Convert to galactic coordinates.
    l = coordinates.galactic.l.degree
    b = coordinates.galactic.b.degree
    theta = (90.-b)*np.pi/180.
    phi = l*np.pi/180.
    print "l, b = %.3f, %.3f" % (l, b)
    print "theta, phi = %.3f, %.3f" % (theta, phi)
    m = load_map(map)

    # from this code
    ebv = healpy.get_interp_val(m, theta, phi)

    # from astroquery
    t = astroquery.get_extinction_table('%.4f %.4f' % (ra, dec))
    if map in [0, 2, 3]:
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
