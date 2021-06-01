"""
Module
------
dataretrieve.py: Data retrieval and output class

Summary
-------
Contains AstroObjectClass which allows collection of lightcurve information for a specific object, based on position.
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

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.io import ascii
from astropy.io.votable import parse
from astroquery.irsa import Irsa
# Set up matplotlib
from matplotlib import pyplot as plt
# from matplotlib import artist
from scipy import stats

from LCExtract import config
from LCExtract.PanSTARRS import ps1cone, getDetections
from LCExtract.coord import CoordClass, to_string
from LCExtract.utilities import Spinner

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
    filterStr = getFilterStr(config.ztf.filters)  # limit filters (requested) to ZTF subset

    status = True
    delim = "%20"
    ra = coordinates.ra_str() + delim
    dec = coordinates.dec_str() + delim
    radius_str = to_string(radius, 5)
    if column_filters is None:
        column_filters = {}

    queryPart = "nph_light_curves"
    pos = "POS=CIRCLE" + delim + ra + dec + radius_str
    bandname = "BANDNAME=" + filterStr
    form = "FORMAT=" + return_type
    badCatFlagsMask = "BAD_CATFLAGS_MASK=32768"

    url_payload = f"{config.ztf.URL}{queryPart}?{pos}&{bandname}&{form}&{badCatFlagsMask}"

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
        return config.badResponse

    memFile = io.BytesIO(siteData.read())

    votable = parse(memFile)
    table = votable.get_first_table().to_table(use_names_over_ids=True)

    if not len(table):  # Check table actually has data in it (i.e. possible no lightcurve data exists)
        return config.badResponse

    tablePD = table.to_pandas()

    fi = pd.Series({"zg": "g", "zr": "r", "zi": "i"})  # map filter ID from ZTF code (used as key in output)
    tablePD['filterID'] = tablePD['filtercode'].map(fi)

    return status, tablePD


def getLightCurveDataPanSTARRS(coords: CoordClass, radius, return_type, column_filters=None):
    """Pan-STARRS light curve data retrieval

    The Pan-STARRs catalog API allows the ability to search the Pan-STARRS catalogs. For additional information
    on the catalogs please visit the Pan-STARRS Data Archive Home Page.

    Ref. https://outerspace.stsci.edu/display/PANSTARRS/Pan-STARRS1+data+archive+home+page


    :param coords: Coordinates of object expressed CoordClass notation in J2000 RA Dec (Decimal) format.
    :type coords: CoordClass
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
        return config.badResponse

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
        return config.badResponse
    else:
        dTab.remove_rows(dTab['psfFlux'] == 0)  # remove any rows where flux is zero
        dTab['psfMag'] = -2.5 * np.log10(dTab['psfFlux']) + 8.90  # convert flux (in Janskys) to magnitude

    return status, dTab.to_pandas()


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

    rt = ('HTML', 'ASCII', 'SVC', 'VOTABLE', 'XML') # not used currently
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
    for key in config.filterSelection:
        if key in statDict.keys():
            if isinstance(statDict[key], (np.float64, np.float32)):
                print(f'{statDict[key]:^{lenVal}.{lenDP}f}', end='')
            elif isinstance(statDict[key], np.int64):
                print(f'{statDict[key]:^{lenVal}}', end='')
        else:
            print(f'{" ":{lenVal}}', end='')
    print()


class AstroObjectClass:
    """Class representing an astronomical object

    """

    def __init__(self, objectName, ra, dec, subtitle=None):
        """
        AstroObjectClass initialises name and object position.
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
        self.pos = CoordClass(ra, dec)

    def preparePlot(self, plotRows):
        fig, ax = plt.subplots(nrows=1, ncols=1, sharex='all', sharey='none')
        ax.set_xlabel('Time [MJD]', fontsize=14)
        ax.set_ylabel('Mag', fontsize=14)
        ax.set_title(self.shortDesc, fontsize=12)

        return fig, ax

    def finalisePlot(self, status, fig, ax):
        if status:
            ymin, ymax = ax.get_ylim()
            if ymax-ymin < 1.0:
                ymid = (ymin + ymax)/2
                ax.set_ylim(ymid-0.5, ymid+0.5)
            ax.set_ylim(reversed(ax.set_ylim()))  # flip the y-axis
            fig.suptitle(f'{self.objectName}', fontsize=16)
            plt.show()


class AODataClass:
    """Class for storing and manipulating the data for an astronomical object"""

    def __init__(self, AO: AstroObjectClass):
        self.table = pd.DataFrame()
        self.samples = {}
        self.mad = {}
        self.SD = {}
        self.median = {}
        self.mean = {}
        self.AO = AO

    def getLightCurveData(self, catalog, radius=None, return_type='VOTABLE'):
        """
        Class method to get light curve data from a particular source

        Radius of search (optional) set first and then catalogue queried through specific method calls. Will be
        overloaded for different catalogs

        :param radius: Cone search radius in arcseconds. Optional
        :type radius: float
        :param catalog: Named tuple representing archive under query
        :type catalog: config.Archive
        :param return_type: Type of return format required. Should be 'VOTABLE'
        :type return_type: str
        :return: Successful extract and data ingested
        :rtype: bool
        """
        if radius is None:
            radiusDeg = config.coneRadius
        else:
            radiusDeg = radius / 3600

        if catalog.name == 'ZTF':
            response = getLightCurveDataZTF(self.AO.pos, radiusDeg, return_type)
        elif catalog.name == 'Pan-STARRS':
            response = getLightCurveDataPanSTARRS(self.AO.pos, radiusDeg, return_type='CSV')
        elif catalog.name == 'PTF':
            response = getLightCurveDataPTF(self.AO.pos, radiusDeg, return_type)
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

    def setMean(self, col_name, group_col):
        """Method to set the median of data

        Value(s) set within the data structure for each individual filter within the data

        :param group_col:
        :type group_col:
        :param col_name: Column name on which to apply the summary, e.g. 'mag'
        :type col_name: str
        """
        self.mean = self.table.groupby(group_col)[col_name].mean()

    def addColourColumn(self, series):
        """Method to add a colour column

        As the data does not have a corresponding colour associated, this is added in a column for
        each row, based on the value of the filter in which the observation was made. This is only used for plots.

        :param series: Column name to use for colour selection
        :type series: str
        """
        c = pd.Series({"g": "green", "r": "red", "i": "indigo", "z": "blue", "y": "black", "R": "orange"})
        self.table['colour'] = self.table[series].map(c)

    def getData(self, archive):
        """Method to encapsulate extraction of data

         Data extracted from catalog and summary statistical analysis carried out.
         All data stored in class structure

        :return: Status of data extract
        :rtype: bool
        """

        if self.getLightCurveData(catalog=archive):
            self.setSamples(archive.magField, archive.filterField)
            self.setMad(archive.magField, archive.filterField)
            self.setSD(archive.magField, archive.filterField)
            self.setMedian(archive.magField, archive.filterField)
            self.setMean(archive.magField, archive.filterField)
            return True
        else:
            return False

    def getTable(self):
        return self.table

    def plot(self, fig, ax, archive):
        """Method to encapsulate the plotting of data

        Sets colour column to distinguish different filters used in data, then sets up plot from given
        X and Y columns. Title is set to object name.

        :param ax:
        :type ax:
        :param x: X-axis title
        :type x: str
        :param y: Y-axis title
        :type y: str
        :param archive: Archive object used for configuration of the plot
        :type archive: config.Archive
        """
        if True:  # TODO Need to sort this out for different catalogs
            self.addColourColumn(archive.filterField)
            colors = self.table['colour']

        # fig, ax = plt.subplots()

        ax.scatter(self.table[archive.timeField], self.table[archive.magField],
                   c=colors, marker=archive.marker)
        # ax.set_ylim(reversed(ax.set_ylim()))  # flip the y-axis
        # plt.xlabel(x, fontsize=14)
        # plt.ylabel(y, fontsize=14)
        # plt.suptitle(f'{self.AO.objectName} {archive.name}', fontsize=16)
        # plt.title(self.AO.shortDesc, fontsize=12)
        # legend1 = ax.legend(*scatter.legend_elements(num=5, c=colors, label=archive.filterField),
        #                     loc="upper right", title="Filters")
        # ax.add_artist(legend1)

        # plt.show()

    def objectOutput(self, archive):
        """Method to encapsulate data output

        Table of summary statistics is sent to console with a plot of data output to plot window.
        """
        print(f"Archive name: {archive.name}")
        print(f'{" ":30}', end='')
        for key in config.filterSelection:
            print(f'{key:^8}', end='')
        print()
        # output filter data line stats to console
        filterLineOut('Samples', self.samples)
        filterLineOut('Median Absolute Deviation', self.mad)
        filterLineOut('Standard Deviation', self.SD)
        filterLineOut('Median', self.median, 2)
        filterLineOut('Mean', self.mean, 2)
        print()
        # plot graph of mag data vs. date
        # self.plot(('mjd', '$mjd$'), ('mag', '$mag$'), 'filtercode')
        # self.plot('$Time [MJD]$', '$mag$', archive)
