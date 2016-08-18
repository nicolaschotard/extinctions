"""
Astronomical info web-based fetchers/parsers.

This code comes from the SNfactory ToolBox (ToolBox.Astro.Fetchers)
"""

from datetime import datetime
import urlparse
import urllib
import re
from astropy.coordinates import SkyCoord


def subdict(dic, keyss=None):
    """Iterate sub-dicts of a dict by key list."""
    if keyss is None:
        keyss = []
    yield dic
    for k in keyss:
        if isinstance(dic, dict) and k in dic:
            dic = dic[k]
            if isinstance(dic, (tuple, list)):
                dic = dic[-1]
            yield dic
        else:
            yield None


xml_keys, xml_info = [], {}


def xml_parser(xml, exclude=None):
    """
    Generalized Expat-XML parser to python dict.

    :param exclude: XML keys to explicitely exclude
    """
    import xml.parsers.expat as expat

    if exclude is None:
        exclude = []
    global xml_keys, xml_info
    xml_keys, xml_info = [], {}

    def start_element(name, attrs):
        xml_keys.append(name)
        if name not in exclude:
            # black magic
            sd = list(subdict(xml_info, xml_keys))
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
        global xml_keys
        if xml_keys[-1] == name:
            xml_keys = xml_keys[:-1]

    def char_data(data):
        if xml_keys[-1] not in exclude and data.strip() != '':
            sd = list(subdict(xml_info, xml_keys))
            sd[-2][xml_keys[-1]] = data.strip()

    p = expat.ParserCreate()
    p.StartElementHandler = start_element
    p.EndElementHandler = end_element
    p.CharacterDataHandler = char_data

    try:
        p.Parse(xml)
    except expat.ExpatError:
        return None

    return xml_info


def sfd_ebmv(ra, dec, equinox=2000, obsepoch=None, service='IRSA'):
    """
    Fetch SFD Milky Way E(B-V) from IRSA or NED.

    `IRSA <http://irsa.ipac.caltech.edu/>`_ (slow)
    `NED <http://ned.ipac.caltech.edu/>`_.

    :param service: can be 'IRSA' or 'NED'
    """
    assert equinox in [1950, 2000], 'equinox should be 1950 or 2000'

    if obsepoch is None:
        obsepoch = datetime.utcnow().year

    ra, dec = SkyCoord(ra, dec, unit='deg').to_string('hmsdms').split()
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
        mwebv = re.search(r'(?<=E\(B-V\)\s=\s)\d+\.\d+(?=\smag)', html)
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
        infos = xml_parser(urllib.urlopen(url).read())
        if infos is not None:
            # the other 'result' are not in mag
            mwebv = infos['results']['result'][0]['statistics']['refPixelValueSFD']
            return float(mwebv.split()[0])
        else:
            return None
