"""
.. _snfactory:

Astronomical info web-based fetchers/parsers
============================================

This code comes from the SNfactory ToolBox (ToolBox.Astro.Fetchers)
"""

from astropy.coordinates import SkyCoord

def subdict(d, keys=None):
    """
    Iterate sub-dicts of a dict by key list
    """
    if keys is None:
        keys = []
    yield d
    for k in keys:
        if isinstance(d, dict) and k in d:
            d = d[k]
            if isinstance(d, (tuple, list)):
                d = d[-1]
            yield d
        else:
            yield None

def xml_parser(xml, exclude=None):
    """
    Generalized Expat-XML parser to python dict

    :param exclude: XML keys to explicitely exclude
    """
    import xml.parsers.expat as expat

    if exclude is None:
        exclude = []
    global keys, info
    keys, info = [], {}

    def start_element(name, attrs):
        global keys, info
        keys.append(name)
        if name not in exclude:
            # black magic
            sd = list(subdict(info, keys))
            # new key
            if sd[-1] is None:
                # add to a sub-dict
                if isinstance(sd[-2], dict):
                    sd[-2].update({name: attrs})
                else:
                    sd[-3].update({name: attrs})
            # already has a dict with the same name at the same level
            elif isinstance(sd[-1], dict):
                if isinstance(sd[-2][name], (tuple, list)):
                    sd[-2][name] += [attrs]
                else:
                    sd[-2][name] = [sd[-2][name], attrs]

    def end_element(name):
        global keys
        if keys[-1] == name:
            keys = keys[:-1]

    def char_data(data):
        global keys, info
        if keys[-1] not in exclude and data.strip() != '':
            sd = list(subdict(info, keys))
            sd[-2][keys[-1]] = data.strip()

    p = expat.ParserCreate()
    p.StartElementHandler = start_element
    p.EndElementHandler = end_element
    p.CharacterDataHandler = char_data

    try:
        p.Parse(xml)
    except expat.ExpatError:
        return None

    return info

def sfd_ebmv(ra, dec, equinox=2000, obsepoch=None, service='IRSA'):
    """
    Fetch SFD Milky Way E(B-V) from `IRSA
    <http://irsa.ipac.caltech.edu/>`_ (slow) or `NED
    <http://ned.ipac.caltech.edu/>`_.

    :param service: can be 'IRSA' or 'NED'
    """
    from datetime import datetime
    import urlparse
    import urllib
    import re

    assert equinox in [1950, 2000], 'equinox should be 1950 or 2000'

    if obsepoch is None:
        obsepoch = datetime.utcnow().year

    c = SkyCoord(ra, dec, unit='deg')
    ra, dec = c.to_string('hmsdms').split()
    equinox = {1950: 'B1950.0', 2000: 'J2000.0'}[equinox]

    if service.lower() == 'ned':

        query = dict(in_csys='Equatorial',
                     in_equinox=equinox,
                     obs_epoch=str(obsepoch),
                     lon=ra.__str__(),
                     lat=dec.__str__(),
                     pa='0.0',
                     out_csys='Equatorial',
                     out_equinox=equinox)

        url = urlparse.urlunparse(('http', 'ned.ipac.caltech.edu',
                                   'cgi-bin/nph-calc', '',
                                   '&'.join(['%s=%s' % (i, j)
                                             for i, j in query.iteritems()]),
                                   ''))
        html = urllib.urlopen(url).read()

        # grep into the HTML
        mwebv = re.search('(?<=E\(B-V\)\s=\s)\d+\.\d+(?=\smag)', html)
        if mwebv is not None:
            return float(mwebv.group(0))
        else:
            return None

    elif service.lower() == 'irsa':

        query = dict(locstr=urllib.quote('%s %s Equatorial %s' %
                                         (ra.__str__(), dec.__str__(), equinox)))

        url = urlparse.urlunparse(('http', 'irsa.ipac.caltech.edu',
                                   'cgi-bin/DUST/nph-dust', '',
                                   '&'.join(['%s=%s' % (i, j)
                                             for i, j in query.iteritems()]),
                                   ''))

        # parse the XML
        info = xml_parser(urllib.urlopen(url).read())
        if info is not None:
            # the other 'result' are not in mag
            mwebv = info['results']['result'][0]['statistics']['refPixelValueSFD']
            return float(mwebv.split()[0])
        else:
            return None
