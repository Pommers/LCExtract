"""
Module
______

Summary
_______

Details
_______

Notes
_____


"""
import json
import sys

import numpy as np
import requests
# imports
# %matplotlib inline - This is needed to insert plots inline in a Jupyter notebook, but not here.
from astropy.io import ascii
from astropy.table import Table

from LCExtract import config

try:  # Python 3.x
    from urllib.parse import quote as urlencode
    from urllib.request import urlretrieve
except ImportError:  # Python 2.x
    from urllib import pathname2url as urlencode
    from urllib import urlretrieve

try:  # Python 3.x
    import http.client as httplib
except ImportError:  # Python 2.x
    import httplib


def ps1cone(ra, dec, radius, table="mean", release="dr1", format="csv", columns=None, verbose=False, **kw):
    """Do a cone search of the PS1 catalog

    **kw: other parameters (e.g., 'nDetections.min':2)
    :param ra: (degrees) J2000 Right Ascension
    :type ra: float
    :param dec: (degrees) J2000 Declination
    :type dec: float
    :param radius: (degrees) Search radius (<= 0.5 degrees)
    :type radius: float
    :param table: mean, stack, or detection
    :type table: str
    :param release: dr1 or dr2
    :type release: str
    :param format: csv, votable, json
    :type format: str
    :param columns: list of column names to include (None means use defaults)
    :type columns: list
    :param baseurl: base URL for the request
    :type baseurl: str
    :param verbose: print info about request
    :type verbose: bool
    :param kw:
    :type kw:
    :return:
    :rtype:
    """

    data = kw.copy()
    data['ra'] = ra
    data['dec'] = dec
    data['radius'] = radius
    return ps1search(table=table, release=release, format=format, columns=columns, verbose=verbose, **data)


def ps1search(table="mean", release="dr1", format="csv", columns=None, verbose=False, **kw):
    """Do a general search of the PS1 catalog (possibly without ra/dec/radius)

    :param table: mean, stack, or detection
    :type table: str
    :param release: dr1 or dr2
    :type release: str
    :param format: csv, votable, json
    :type format: str
    :param columns: list of column names to include (None means use defaults)
    :type columns: list
    :param verbose: print info about request
    :type verbose: bool
    :param kw: other parameters (e.g., 'nDetections.min':2).  Note this is required!
    :type kw: str
    :return: search response from MAST based Pan-STARRS archive
    :rtype: str
    """

    baseurl = config.panstarrs.URL  # base URL for the Pan-STARRS request from MAST in config file
    data = kw.copy()
    if not data:
        raise ValueError("You must specify some parameters for search")
    checklegal(table, release)
    if format not in ("csv", "votable", "json"):
        raise ValueError("Bad value for format")
    url = f"{baseurl}/{release}/{table}.{format}"
    url = "{baseurl}/{release}/{table}.{format}".format(**locals())
    if columns:
        # check that column values are legal
        # create a dictionary to speed this up
        dcols = {}
        for col in ps1metadata(table, release)['name']:
            dcols[col.lower()] = 1
        badcols = []
        for col in columns:
            if col.lower().strip() not in dcols:
                badcols.append(col)
        if badcols:
            raise ValueError('Some columns not found in table: {}'.format(', '.join(badcols)))
        # two different ways to specify a list of column values in the API
        # data['columns'] = columns
        data['columns'] = '[{}]'.format(','.join(columns))

    # either get or post works
    #    r = requests.post(url, data=data)
    r = requests.get(url, params=data)

    if verbose:
        print(r.url)
    r.raise_for_status()
    if format == "json":
        return r.json()
    else:
        return r.text


def checklegal(table, release):
    """Checks if this combination of table and release is acceptable

    Raises a VelueError exception if there is problem
    """

    releaselist = ("dr1", "dr2")
    if release not in ("dr1", "dr2"):
        raise ValueError("Bad value for release (must be one of {})".format(', '.join(releaselist)))
    if release == "dr1":
        tablelist = ("mean", "stack")
    else:
        tablelist = ("mean", "stack", "detection")
    if table not in tablelist:
        raise ValueError("Bad value for table (for {} must be one of {})".format(release, ", ".join(tablelist)))


def ps1metadata(table="mean", release="dr1",
                baseurl="https://catalogs.mast.stsci.edu/api/v0.1/panstarrs"):
    """Return metadata for the specified catalog and table

    :param table: mean, stack, or detection
    :type table: str
    :param release: dr1 or dr2
    :type release: str
    :param baseurl: base URL for the request
    :type baseurl: str

    :return Returns an astropy table with columns name, type, description
    :rtype astropy table

    """

    checklegal(table, release)
    url = "{baseurl}/{release}/{table}/metadata".format(**locals())
    r = requests.get(url)
    r.raise_for_status()
    v = r.json()
    # convert to astropy table
    tab = Table(rows=[(x['name'], x['type'], x['description']) for x in v],
                names=('name', 'type', 'description'))
    return tab


def mastQuery(request):
    """Perform a MAST query.

    Parameters
    ----------
    request (dictionary): The MAST request json object

    Returns head,content where head is the response HTTP headers, and content is the returned data
    """

    server = 'mast.stsci.edu'

    # Grab Python Version
    version = ".".join(map(str, sys.version_info[:3]))

    # Create Http Header Variables
    headers = {"Content-type": "application/x-www-form-urlencoded",
               "Accept": "text/plain",
               "User-agent": "python-requests/" + version}

    # Encoding the request as a json string
    requestString = json.dumps(request)
    requestString = urlencode(requestString)

    # opening the https connection
    conn = httplib.HTTPSConnection(server)

    # Making the query
    conn.request("POST", "/api/v0/invoke", "request=" + requestString, headers)

    # Getting the response
    resp = conn.getresponse()
    head = resp.getheaders()
    content = resp.read().decode('utf-8')

    # Close the https connection
    conn.close()

    return head, content


def resolve(name):
    """Get the RA and Dec for an object using the MAST name resolver

    :param name: Name of object
    :type name: str
    :return: RA, Dec tuple with position
    :rtype: tuple
    """

    resolverRequest = {'service': 'Mast.Name.Lookup',
                       'params': {'input': name,
                                  'format': 'json'
                                  },
                       }
    headers, resolvedObjectString = mastQuery(resolverRequest)
    resolvedObject = json.loads(resolvedObjectString)
    # The resolver returns a variety of information about the resolved object,
    # however for our purposes all we need are the RA and Dec
    try:
        objRa = resolvedObject['resolvedCoordinate'][0]['ra']
        objDec = resolvedObject['resolvedCoordinate'][0]['decl']
    except IndexError as e:
        raise ValueError("Unknown object '{}'".format(name))
    return objRa, objDec


def addFilter(dTab: Table):
    """Add filter name as column in detection table by translating filterID

    This modifies the table in place.  If the 'filter' column already exists,
    the table is returned unchanged.
    """
    if 'filtercode' not in dTab.colnames:
        # the filterID value goes from 1 to 5 for grizy
        id2filter = np.array(list('grizy'))
        dTab['filtercode'] = id2filter[(dTab['filterID'] - 1).data]
        # TODO need to remove filter data which is not requested
    return dTab


def getDetections(tab):
    """Extract objects from the Detection table

    Extract information for object (via ID) with the same object ID from the Detection table, which contains all
    the individual measurements for a source

    :param tab: Table of sources about a position
    :type tab: Table
    :return: Detections for a single source from the detections table
    :rtype: Table
    """
    objID = tab['objID'][0]  # only first object selected...
    filters = [i + 1 for i, n in enumerate('grizy') if n in config.filterSelection]
    dConstraints = {'objID': objID, 'filterID': filters}
    dColumns = ("""objID,detectID,filterID,obsTime,ra,dec,psfFlux,psfFluxErr,psfMajorFWHM,psfMinorFWHM,
                psfQfPerfect,apFlux,apFluxErr,infoFlag,infoFlag2,infoFlag3""").split(',')
    # strip blanks and weed out blank and commented-out values
    dColumns = [x.strip() for x in dColumns]
    dColumns = [x for x in dColumns if x and not x.startswith('#')]

    dResults = ps1search(table='detection', release='dr2', columns=dColumns, **dConstraints)
    dTab = addFilter(ascii.read(dResults))
    dTab.sort('obsTime')
    return dTab
