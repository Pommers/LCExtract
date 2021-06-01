from urllib.error import HTTPError

import pandas as pd
from astropy import units as u
from astroquery.irsa import Irsa

from LCExtract import config
from LCExtract.coord import CoordClass
from LCExtract.filter import getFilterStr
from LCExtract.utilities import Spinner


def getLightCurveDataPTF(coordinates: CoordClass, radius,
                         return_type, column_filters=None):
    """Palomar Transient factory light curve data retrieval

    IRSA provides access to the PTF collection of lightcurve data through an application program interface (API).
    Search, restriction, and formatting parameters are specified in an HTTP URL. The output is a table in the
    requested format containing lightcurve data satisfying the search constraints.

    Ref. https://irsa.ipac.caltech.edu/applications/Gator/GatorAid/irsa/catsearch.html

    :param coordinates: Coordinates of object expressed CoordClass notation in J2000 RA Dec (Decimal) format.
    :type coordinates: CoordClass
    :param radius: Radius of cone search ** in degrees ** for passing to PTF
    :type radius: float
    :param return_type: For selection of different return types, e.g. "VOTABLE" (Default), "HTML", "CSV"
    :type return_type: str
    :param column_filters: Not used currently
    :returns:
        (boolean) Valid data return
        (DataFrame) Data payload
    :rtype: tuple

    """

    filterStr = getFilterStr(config.ptf.filters)  # limit filters (requested) to PTF subset

    status = True

    rt = ('HTML', 'ASCII', 'SVC', 'VOTABLE', 'XML')  # not used currently
    if column_filters is None:
        column_filters = {}

    print('Requesting data from Palomar Transient Factory. Please wait ... ', end='')
    with Spinner():
        try:
            votable = Irsa.query_region(f"{coordinates.getRA()}, {coordinates.getDEC()}", catalog="ptf_lightcurves",
                                        spatial="Cone", radius=radius * u.deg, verbose=False)
            print(f'\r{" ":65}\r ', end='')
        except HTTPError as err:
            if err.code == 400:
                print('Sorry. Could not complete request.')
            else:
                raise

    if not len(votable):  # Check table actually has data in it (i.e. possible no lightcurve data exists)
        return config.badResponse

    tablePD = votable.to_pandas()

    # Filter constraint - fid=1 (g filter) or fid=2 (R filter)
    fi = pd.Series({1: "g", 2: "R"})  # map filter ID from ZTF code (used as key in output)
    tablePD['filterID'] = tablePD['fid'].map(fi)
    tablePD = tablePD.loc[tablePD['filterID'].isin(list(config.filterSelection))]

    return status, tablePD