"""
From http://argonaut.skymaps.info/usage#function-call.

With a few modifications to handle large query.
"""

import json
import requests


def query(lon, lat, coordsys='gal', mode='full', limit=500000):
    """
    Send a line-of-sight reddening query to the Argonaut web server.

    lon, lat: longitude and latitude, in degrees.
    coordsys: 'gal' for Galactic, 'equ' for Equatorial (J2000).
    mode: 'full', 'lite' or 'sfd'

    In 'full' mode, outputs a dictionary containing, among other things:

    - 'distmod':    The distance moduli that define the distance bins.
    - 'best':       The best-fit (maximum proability density)
                  line-of-sight reddening, in units of SFD-equivalent
                  E(B-V), to each distance modulus in 'distmod.' See
                  Schlafly & Finkbeiner (2011) for a definition of the
                  reddening vector (use R_V = 3.1).
    - 'samples':    Samples of the line-of-sight reddening, drawn from
                  the probability density on reddening profiles.
    - 'success':    1 if the query succeeded, and 0 otherwise.
    - 'converged':  1 if the line-of-sight reddening fit converged, and
                  0 otherwise.
    - 'n_stars':    # of stars used to fit the line-of-sight reddening.
    - 'DM_reliable_min':  Minimum reliable distance modulus in pixel.
    - 'DM_reliable_max':  Maximum reliable distance modulus in pixel.

    Less information is returned in 'lite' mode, while in 'sfd' mode,
    the Schlegel, Finkbeiner & Davis (1998) E(B-V) is returned.
    """
    # make sure we have list
    if type(lon) == float:
        lon, lat = [lon], [lat]

    # Make sure to have less than 500000 objects (the limit).
    # Cut the list in smaller pieces if that is the case.
    if len(lon) >= limit:
        dicts = [query(loni, lati, coordsys=coordsys, mode=mode)
                 for loni, lati in zip(chunk(lon, limit - 1),
                                       chunk(lat, limit - 1))]
        for dic in dicts[1:]:
            for k in dic:
                dicts[0][k].extend(dic[k])
        return dicts[0]

    if coordsys.lower() in ['gal', 'g']:
        payload = {'mode': mode, 'l': lon, 'b': lat}
    elif coordsys.lower() in ['equ', 'e']:
        payload = {'mode': mode, 'ra': lon, 'dec': lat}
    else:
        raise ValueError("coordsys '{0}' not understood.".format(coordsys))

    req = requests.post('http://argonaut.skymaps.info/gal-lb-query-light',
                        data=json.dumps(payload),
                        headers={'content-type': 'application/json'})

    try:
        req.raise_for_status()
    except requests.exceptions.HTTPError as excep:
        print 'Response received from Argonaut:'
        print req.text
        raise excep

    return json.loads(req.text)


def chunk(ilist, length):
    """Divide a list into smaller lists of maximal length."""
    return [ilist[i:i + length] for i in range(0, len(ilist), length)]
