"""
Module
------
dataretrieve.py: Data retrieval and output class

Summary
-------
Contains DataClass which allows collection of lightcurve information for a specific object, based on position.
Implements specific methods for download of information from archives

1. Zwicky Transient Facility
2. Pan-STARRS


Notes
-----

"""
import io
import re
from urllib.error import HTTPError
from urllib.request import urlopen

import numpy as np
import pandas as pd
from astropy.io.votable import parse
from astropy.io import ascii
# Set up matplotlib
from matplotlib import pyplot as plt
from scipy import stats

from LCExtract import config
from LCExtract.coord import CoordClass, to_string
from LCExtract.utilities import Spinner
from LCExtract.PanSTARRS import ps1cone, ps1search, ps1metadata
from LCExtract.PanSTARRS import checklegal, mastQuery, resolve, getDetections

"""filter dict and list for reference in output iteration"""
ZTFfilters = {"zg": 0, "zr": 1, "zi": 2}
filterKey = list(ZTFfilters)


# TODO Need to implement a way to consolidate filters from different sources


def getFilterStr(avail: str, delim=','):
    """Return a subset of filters requested as appropriate to archive

    :param avail: available filters for archive facility (e.g. 'gri')
    :type avail: str
    :param delim: response delimiter (default ',')
    :type delim: str
    :return: filter subset based on request (e.g. 'g,r')
    :rtype: str
    """
    temp = avail if not config.filterSelection else re.findall('[' + config.filterSelection + ']', avail)
    return delim.join(temp)


def getLightCurveDataZTF(coordinates: CoordClass, radius,
                         return_type, column_filters=None):
    """Zwicky Transient facility light curve data retrieval

    IRSA provides access to the ZTF collection of lightcurve data through an application program interface (API).
    Search, restriction, and formatting parameters are specified in an HTTP URL. The output is a table in the
    requested format containing lightcurve data satisfying the search constraints.

    Ref. https://irsa.ipac.caltech.edu/docs/program_interface/ztf_lightcurve_api.html

    :param coordinates: Coordinates of object expressed CoordClass notation in J2000 RA Dec (Decimal) format.
    :type coordinates: CoordClass
    :param radius: Radius of cone search ** in degrees ** for passing to ZTF
    :type radius: float
    :param return_type: For selection of different return types, e.g. "VOTABLE" (Default), "HTML", "CSV"
    :type return_type: str
    :param column_filters: Not used currently
    :returns:
        (boolean) Valid data return
        (DataFrame) Data payload
    :rtype: tuple

    """
    filterStr = getFilterStr('gri')  # limit filters (requested) to ZTF subset

    status = True
    delim = "%20"
    ra = coordinates.ra_str() + delim
    dec = coordinates.dec_str() + delim
    radius_str = to_string(radius, 4)
    if column_filters is None:
        column_filters = {}

    queryPart = "nph_light_curves"
    pos = "POS=CIRCLE" + delim + ra + dec + radius_str
    bandname = "BANDNAME=" + filterStr
    format = "FORMAT=" + return_type
    badCatFlagsMask = "BAD_CATFLAGS_MASK=32768"

    url_payload = f"{config.baseURL.ZTF}{queryPart}?{pos}&{bandname}&{format}&{badCatFlagsMask}"

    # establish http connection
    # http = urllib3.PoolManager()
    # siteData = http.request('GET', url_payload)
    print('Requesting data from Zwicky Transient Facility. Please wait ... ', end='')
    with Spinner():
        try:
            siteData = urlopen(url_payload)
            print(f'\r{" ":66}\r ', end='')
            # print(' ', end='')
        except HTTPError as err:
            if err.code == 400:
                print('Sorry. Could not complete request.')
            else:
                raise

    if siteData.status != 200:  # Ensure good response is received back from IRSA
        status = False

    memFile = io.BytesIO(siteData.read())

    votable = parse(memFile)
    table = votable.get_first_table().to_table(use_names_over_ids=True)

    if not len(table):  # Check table actually has data in it (i.e. possible no lightcurve data exists)
        status = False

    return status, table.to_pandas()


def getLightCurveDataPanSTARRS(coords: CoordClass, radius, return_type, column_filters=None):
    """Pan-STARRS light curve data retrieval

    The Pan-STARRs catalog API allows the ability to search the Pan-STARRS catalogs. For additional information
    on the catalogs please visit the Pan-STARRS Data Archive Home Page.

    Ref. https://outerspace.stsci.edu/display/PANSTARRS/Pan-STARRS1+data+archive+home+page


    :param coordinates: Coordinates of object expressed CoordClass notation in J2000 RA Dec (Decimal) format.
    :type coordinates: CoordClass
    :param radius: Radius of cone search ** in degrees ** for passing to Pan-STARRS
    :type radius: float
    :param return_type: For selection of different return types, e.g. "VOTABLE" (Default), "HTML", "CSV"
    :type return_type: str
    :param column_filters: Not used currently
    :returns:
        (boolean) Valid data return
        (DataFrame) Data payload
    :rtype: tuple

    """

    constraints = {'nDetections.gt': 1}
    # set columns to return by default
    # strip blanks and weed out blank and commented-out values
    columns = """objID,raMean,decMean,nDetections,ng,nr,ni,nz,ny,
        gMeanPSFMag,rMeanPSFMag,iMeanPSFMag,zMeanPSFMag,yMeanPSFMag""".split(',')
    columns = [x.strip() for x in columns]
    columns = [x for x in columns if x and not x.startswith('#')]

    # limit filters (requested) to PanSTARRS subset
    filterStr = getFilterStr('grizy')

    status = True
    if column_filters is None:
        column_filters = {}

    print('Searching for object in Pan-STARRS archive (MAST). Please wait ... ', end='')
    with Spinner():
        try:
            # perform a cone search about coordinates to get detections
            results = ps1cone(coords.getRA(), coords.getDEC(), radius, release='dr2', columns=columns, **constraints)
            print(f'\r{" ":68}\r ', end='')
        except HTTPError as err:
            if err.code == 400:
                print('Sorry. Could not complete request.')
            else:
                raise

    if not results:
        return False

    # convert to table
    tab = ascii.read(results)

    # improve the format
    for filter in 'grizy':
        col = filter + 'MeanPSFMag'
        tab[col].format = ".4f"  # (only for printing?)
        tab[col][tab[col] == -999.0] = np.nan  # set to nan if -999 before analysis

    print('Searching for object detections. Please wait ... ', end='')
    with Spinner():
        try:
            # get individual detections for first object in the list
            dTab = getDetections(tab)
            print(f'\r{" ":50}\r ', end='')

        except HTTPError as err:
            if err.code == 400:
                print('Sorry. Could not complete request.')
            else:
                raise

    if not len(dTab):  # Check table actually has data in it (i.e. possible no lightcurve data exists)
        status = False
    else:
        dTab['mag'] = -2.5 * np.log10(dTab['psfFlux']) + 8.90

    return status, dTab.to_pandas()


def filterLineOut(statStr, statDict, lenDP=3, lenStr=30, lenVal=8, valueType=float):
    """Output line of individual filter data to the console

    e.g. "Median Absolute Deviation      0.039   0.026   0.024  "

    :param statStr: String describing filter output
    :type statStr: str
    :param statDict: Dictionary of filter / summary statistic pairs
    :type statDict: dict
    :param lenDP: Number of decimal places for the value display (Optional, Default=3)
    :type lenDP: int
    :param lenStr: Length of the stat summary string (Optional, Default=30)
    :type lenStr: int
    :param lenVal: Total length of the value display (Optional, Default=8)
    :type lenVal: int
    """
    print(f'{statStr:{lenStr}}', end='')
    for key in config.filterSelection:
        if key in statDict.keys():
            if isinstance(statDict[key], np.float64):
                print(f'{statDict[key]:^{lenVal}.{lenDP}f}', end='')
            elif isinstance(statDict[key], np.int64):
                print(f'{statDict[key]:^{lenVal}}', end='')
        else:
            print(f'{" ":{lenVal}}', end='')
    print()


class DataClass:
    """Class representing the data for an astronomical object

    """

    def __init__(self, objectName, ra, dec, subtitle=None):
        """
        DataClass initialises name and object position.
        Other default parameters for searches also set as well as structure for data.

        :param subtitle: A string used for information about the object, e.g. location, type/min-desc,
        effective radius, distance. Currently used for subtitle of graph page.
        :type subtitle:
        :param objectName: The name of the object of interest, used to describe it in output. Not used for comparison.
        :type objectName: str
        :param ra: RA value in degrees for the object position
        :type ra: float
        :param dec: DEC value in degrees for the object position
        :type dec: float
        """

        self.objectName = objectName
        self.shortDesc = subtitle
        self.radius = 1 / 3600  # 1 arcseconds
        self.table = pd.DataFrame()
        self.data = np
        self.pos = CoordClass(ra, dec)
        self.samples = {}
        self.mad = {}
        self.SD = {}
        self.median = {}
        self.filters = 'g,r,i,z'

    def getLightCurveData(self, catalog, radius=None, return_type='VOTABLE'):
        """
        Class method to get light curve data from a particular source

        Radius of search (optional) set first and then catalogue queried through specific method calls. Will be
        overloaded for different catalogs

        :param radius: Cone search radius in arcseconds. Optional
        :type radius: float
        :param catalog:
        :type catalog: str
        :param return_type: Type of return format required. Should be 'VOTABLE'
        :type return_type: str
        :return: Successful extract and data ingested
        :rtype: bool
        """
        if radius is None:
            radiusDeg = self.radius
        else:
            radiusDeg = radius / 3600

        if catalog[0] == 'ZTF':
            response = getLightCurveDataZTF(self.pos, radiusDeg, return_type)
        elif catalog[0] == 'PanSTARRS':
            response = getLightCurveDataPanSTARRS(self.pos, radiusDeg, return_type='CSV')
        else:
            return False

        if response[0]:
            self.table = response[1]
            return True
        else:
            return False

    def getCol(self, col_name):
        return self.table[col_name]

    def setSamples(self, col_name, group_col):
        self.samples = self.table.groupby(group_col)[col_name].count()

    def setMad(self, col_name, group_col):
        """Method to set the median absolute deviation

        Value(s) set within the data structure for each individual filter within the data

        :param group_col: column name in series on which to group for filter data
        :type group_col: str
        :param col_name: Column name on which to apply the summary, e.g. 'mag'
        :type col_name: str
        """
        series = self.table.groupby(group_col)[col_name]
        for name, group in series:
            self.mad[name] = stats.median_abs_deviation(series.get_group(name))

        # TODO Need to make grouping more generic, i.e. if 'filtercode' is not the column name, or if only one filter
        #  value exists for a data set. This applies to all summary statistics below.

    def setSD(self, col_name, group_col):
        """Method to set the standard deviation

        Value(s) set within the data structure for each individual filter within the data

        :param group_col:
        :type group_col:
        :param col_name: Column name on which to apply the summary, e.g. 'mag'
        :type col_name: str
        """
        self.SD = self.table.groupby(group_col)[col_name].std()

    def setMedian(self, col_name, group_col):
        """Method to set the median of data

        Value(s) set within the data structure for each individual filter within the data

        :param group_col:
        :type group_col:
        :param col_name: Column name on which to apply the summary, e.g. 'mag'
        :type col_name: str
        """
        self.median = self.table.groupby(group_col)[col_name].median()

    def addColourColumn(self, series):
        """Method to add a colour column

        As the data does not have a corresponding colour associated, this is added in a column for
        each row, based on the value of the filter in which the observation was made. This is only used for plots.

        :param series: Column name to use for colour selection
        :type series: str
        """
        c = pd.Series({"g": "green", "r": "red", "i": "indigo", "z": "blue", "y": "black"})
        self.table['colour'] = self.table[series].map(c)

    def plot(self, x, y, series):
        """Method to encapsulate the plotting of data

        Sets colour column to distinguish different filters used in data, then sets up plot from given
        X and Y columns. Title is set to object name.

        :param x: X-axis column name and title
        :type x: tuple
        :param y: Y-axis column name and title
        :type y: tuple
        :param series: Column name to use to set colour of series
        :type series: str

        """
        if True:  # TODO Need to sort this out for different catalogs
            self.addColourColumn(series)
            colors = self.table['colour']

        plt.scatter(self.table[x[0]], self.table[y[0]], c=colors, marker='*')
        plt.ylim(reversed(plt.ylim()))  # flip the y-axis
        plt.xlabel(x[1], fontsize=14)
        plt.ylabel(y[1], fontsize=14)
        plt.suptitle(self.objectName, fontsize=16)
        plt.title(self.shortDesc, fontsize=12)
        plt.show()

    def getData(self, archive):
        """Method to encapsulate extraction of data

         Data extracted from catalog and summary statistical analysis carried out.
         All data stored in class structure

        :return: Status of data extract
        :rtype: bool
        """

        status = True
        if self.getLightCurveData(catalog=archive):
            self.setSamples('mag', 'filtercode')
            self.setMad('mag', 'filtercode')
            self.setSD('mag', 'filtercode')
            self.setMedian('mag', 'filtercode')
            return status
        else:
            return False

    def objectOutput(self, archive):
        """Method to encapsulate data output

        Table of summary statistics is sent to console with a plot of data output to plot window.
        """
        print(f"Object name: {self.objectName} - summary statistics")
        print(f'{" ":30}', end='')
        for key in config.filterSelection:
            print(f'{key:^8}', end='')
        print()
        # output filter data line stats to console
        filterLineOut('Samples', self.samples)
        filterLineOut('Median Absolute Deviation', self.mad)
        filterLineOut('Standard Deviation', self.SD)
        filterLineOut('Median', self.median, 2)
        # plot graph of mag data vs. date
        # self.plot(('mjd', '$mjd$'), ('mag', '$mag$'), 'filtercode')
        self.plot((archive.timeField, '$Time [MJD]$'), (archive.magField, '$mag$'), 'filtercode')
