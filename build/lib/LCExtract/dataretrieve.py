"""
Module
------
dataretrieve.py: Data retrieval and output class

Summary
-------
Contains DataClass which allows collection of lightcurve information for a specific object, based on position.
Implements specific methods for download of information from archives

1. Zwicky Transient Facility


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
# Set up matplotlib
from matplotlib import pyplot as plt
from scipy import stats

from LCExtract import config
from LCExtract.coord import CoordClass, to_string
from LCExtract.utilities import Spinner

"""filter dict and list for reference in output iteration"""
filters = {"zg": 0, "zr": 1, "zi": 2}
filterKey = list(filters)


# TODO Need to implement a way to consolidate filters from different sources


def getFilterStr(avail: str):
    temp = avail if not config.filterSelection else re.findall('[' + config.filterSelection + ']', avail)
    c = len(temp)
    filt = ''
    for char in range(c):
        filt += temp[char]  # transform to comma separation
        if char != c - 1:
            filt += ','
    return filt


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

    base_url = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/"
    query_url_part = "nph_light_curves?"
    url_pos = "POS=CIRCLE" + delim + ra + dec + radius_str
    url_bandname = "&BANDNAME=" + filterStr
    url_format = "&FORMAT=" + return_type
    url_badCatFlagsMask = "&BAD_CATFLAGS_MASK=32768"
    url_payload = base_url + query_url_part + \
                  url_pos + \
                  url_bandname + \
                  url_format + \
                  url_badCatFlagsMask

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


def getLightCurveDataPanSTARRS(coordinates: CoordClass, radius, return_type, column_filters=None):
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
    # TODO Note yet completed. Access to Pan-STARRS data requires more development

    status = True
    delim = "&"
    ra = coordinates.ra_str() + delim
    dec = coordinates.dec_str() + delim
    radius_str = to_string(radius, 4) + delim
    if column_filters is None:
        column_filters = {}

    base_url = "https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/"
    query_url_part = "dr1/mean?"
    url_pos = ra + dec + radius_str
    url_bandname = "&BANDNAME=g,r,i"
    url_format = "&FORMAT=" + return_type
    url_badCatFlagsMask = "&BAD_CATFLAGS_MASK=32768"
    url_payload = base_url + query_url_part + \
                  url_pos + \
                  url_bandname + \
                  url_format + \
                  url_badCatFlagsMask

    # establish http connection
    # http = urllib3.PoolManager()
    # siteData = http.request('GET', url_payload)
    siteData = urlopen(url_payload)
    if siteData.status != 200:
        status = False

    memFile = io.BytesIO(siteData.read())

    votable = parse(memFile)
    table = votable.get_first_table().to_table(use_names_over_ids=True)

    return status, table.to_pandas()


def filterLineOut(statStr, statDict, lenDP=3, lenStr=30, lenVal=8):
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
    for key in filters:
        if key in statDict.keys():
            print(f'{statDict[key]:^{lenVal}.{lenDP}f}', end='')
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
        self.mad = {}
        self.SD = {}
        self.median = {}
        self.filters = 'g,r,i,z'

    def getLightCurveData(self, radius=None, catalog='ZTF', return_type='VOTABLE'):
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

        if catalog == 'ZTF':
            response = getLightCurveDataZTF(self.pos, radiusDeg, return_type)
            if response[0]:
                self.table = response[1]
                return True
        else:
            return False

    def getCol(self, col_name):
        return self.table[col_name]

    def setMad(self, col_name):
        """Method to set the median absolute deviation

        Value(s) set within the data structure for each individual filter within the data

        :param col_name: Column name on which to apply the summary, e.g. 'mag'
        :type col_name: str
        """
        series = self.table.groupby("filtercode")[col_name]
        for name, group in series:
            self.mad[name] = stats.median_abs_deviation(series.get_group(name))

        # TODO Need to make grouping more generic, i.e. if 'filtercode' is not the column name, or if only one filter
        #  value exists for a data set. This applies to all summary statistics below.

    def setSD(self, col_name):
        """Method to set the standard deviation

        Value(s) set within the data structure for each individual filter within the data

        :param col_name: Column name on which to apply the summary, e.g. 'mag'
        :type col_name: str
        """
        self.SD = self.table.groupby("filtercode")[col_name].std()

    def setMedian(self, col_name):
        """Method to set the median of data

        Value(s) set within the data structure for each individual filter within the data

        :param col_name: Column name on which to apply the summary, e.g. 'mag'
        :type col_name: str
        """
        self.median = self.table.groupby("filtercode")[col_name].median()

    def addColourColumn(self, series):
        """Method to add a colour column

        As the data does not have a corresponding colour associated, this is added in a column for
        each row, based on the value of the filter in which the observation was made. This is only used for plots.

        :param series: Column name to use for colour selection
        :type series: str
        """
        c = pd.Series({"zg": "green", "zr": "red", "zi": "indigo"})
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

        plt.scatter(self.table[x[0]], self.table[y[0]], c=colors)
        plt.ylim(reversed(plt.ylim()))  # flip the y-axis
        plt.xlabel(x[1], fontsize=14)
        plt.ylabel(y[1], fontsize=14)
        plt.suptitle(self.objectName, fontsize=16)
        plt.title(self.shortDesc, fontsize=12)
        plt.show()

    def getData(self):
        """Method to encapsulate extraction of data

         Data extracted from catalog and summary statistical analysis carried out.
         All data stored in class structure

        :return: Status of data extract
        :rtype: bool
        """

        status = True
        if self.getLightCurveData():
            self.setMad('mag')
            self.setSD('mag')
            self.setMedian('mag')
            return status
        else:
            return False

    def objectOutput(self):
        """Method to encapsulate data output

        Table of summary statistics is sent to console with a plot of data output to plot window.
        """
        print(f"Object name: {self.objectName} - summary statistics")
        print(f'{" ":30}{filterKey[0]:^8}{filterKey[1]:^8}{filterKey[2]:^8}')
        # output filter data line stats to console
        filterLineOut('Median Absolute Deviation', self.mad)
        filterLineOut('Standard Deviation', self.SD)
        filterLineOut('Median', self.median, 2)
        # plot graph of mag data vs. date
        self.plot(('mjd', '$mjd$'), ('mag', '$mag$'), 'filtercode')
