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
3. Palomar Transient Factory

Notes
-----

"""
import datetime
import os
import re

import numpy as np
import pandas as pd
from astropy.table import Table

# Set up matplotlib
from matplotlib import pyplot as plt
# from matplotlib import artist
from scipy import stats
# from johansen import Johansen
from astropy.wcs import WCS
from astropy.io import fits
from astropy import time as astroTime

from LCExtract import config
from LCExtract.config import LClog
from LCExtract.coord import CoordClass
from LCExtract.filter import filterLineOut
from LCExtract.utilities import threeSigma
from LCExtract.ztf import getLightCurveDataZTF, getOIDZTFinfo, OIDListFile, getZTFImage, getZTFRefImage
from LCExtract.ptf import getLightCurveDataPTF
from LCExtract.PanSTARRS import getLightCurveDataPanSTARRS


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
        self.plotRows = 0
        self.filtersToPlot = ''

    def incrPlotRows(self):
        self.plotRows += 1

    def plotRowsTrue(self):
        return self.plotRows > 0

    def append_filtersToPlot(self, filtersInArchive):
        self.filtersToPlot += ''.join(filtersInArchive)

    def get_filtersInAO(self):
        if self.filtersToPlot != '':
            return ''.join(re.findall(f'[{self.filtersToPlot}]', config.filterSelection))
        else:
            return ''

    def preparePlot(self, filters):
        fig, ax = plt.subplots(num=1, nrows=len(filters), ncols=1, sharex='all', sharey='none', figsize=(8, 11),
                               dpi=300)
        fig.suptitle(f'{self.objectName}', fontsize=16)

        return fig, ax

    def finalisePlot(self, fig, ax, filters):
        if len(filters) > 1:
            ax[len(filters) - 1].set_xlabel('Date', fontsize=14)
            ax[0].xaxis.set_major_locator(plt.MaxNLocator(6))

            ax[0].set_title(f'{self.shortDesc}', fontsize=12)
            for c in ax:
                r = c.get_subplotspec().rowspan.start
                c.set_ylabel(f'{filters[r]} [Mag]', fontsize=12)
                ymin, ymax = c.get_ylim()
                if ymax - ymin < 1.0:
                    ymid = (ymin + ymax) / 2
                    c.set_ylim(ymid - 0.5, ymid + 0.5)
                c.set_ylim(reversed(c.set_ylim()))  # flip the y-axis
                c.legend()
        else:
            ax.set_xlabel('Date', fontsize=12)
            ax.xaxis.set_major_locator(plt.MaxNLocator(6))

            ax.set_ylabel(f'{filters}[Mag]', fontsize=12)
            ax.set_title(f'{self.shortDesc}', fontsize=12)
            ymin, ymax = ax.get_ylim()
            if ymax - ymin < 1.0:
                ymid = (ymin + ymax) / 2
                ax.set_ylim(ymid - 0.5, ymid + 0.5)
            ax.set_ylim(reversed(ax.set_ylim()))  # flip the y-axis
            ax.legend()

        plotName = self.objectName if self.objectName != '' else datetime.datetime.now().strftime("%y%m%d%H%M%S")
        # plot has been created - now save (as PNG) and output to screen (if required)
        plt.savefig(f'data/plots/{plotName}_{filters}.png')
        if config.plotToScreen:
            plt.show()
        plt.close(fig=1)


def addColumns(r, filters):
    r.add_column(False, name='ZTFObject')
    r.add_column(False, name='ZTFOutliers')
    subSetFilters = ''.join(sorted(set(config.ztf.filters) & set(filters), key=config.ztf.filters.index))
    for f in subSetFilters:
        r.add_column(0, name=f'ZTF{f}oid')
        r.add_column('', name=f'ZTF{f}filterID')
        r.add_column(0, name=f'ZTF{f}Samples')
        r.add_column(0.0, name=f'ZTF{f}MAD')
        r.add_column(0.0, name=f'ZTF{f}SD')
        r.add_column(0.0, name=f'ZTF{f}Median')
        r.add_column(0.0, name=f'ZTF{f}Mean')
        r.add_column(False, name=f'Stat{f}Included')
        r.add_column(0.0, name=f'ZTF{f}RefMag')
        r.add_column(0.0, name=f'ZTF{f}med_16pc_diff')
        r.add_column(0.0, name=f'ZTF{f}deltaTMin')
        r.add_column(0.0, name=f'ZTF{f}deltaTMean')
        r.add_column(0.0, name=f'ZTF{f}deltaTMedian')
        r.add_column(0.0, name=f'ZTF{f}deltaTMax')
        r.add_column(0.0, name=f'ZTF{f}deltaTSD')


class AODataClass:
    """Class for storing and manipulating the data for an astronomical object"""

    def __init__(self, AO: AstroObjectClass, archive: config.Archive):
        self.AO = AO
        self.archive = archive
        self.table = pd.DataFrame()  # This dataframe contains the response from the archive server
        self.samples = {}  # The count of samples returned for each filter
        self.mad = {}  # The median absolute deviation of data from each filter
        self.SD = {}  # The std Dev of data from each filter
        self.median = {}  # the median of data from each filter
        self.mean = {}  # the mean of data from each filter
        self.filtersReturned = {}  # df of filter id based on data returned
        self.id = {}
        self.aRefmag = {}
        self.altSD = {}
        self.sdssRefmag = {}
        self.filename = ''
        self.dataStatus = False
        self.sigma3 = {'g': 0.0, 'r': 0.0}
        self.timeDiff = {}
        self.timeMin = {}
        self.timeMax = {}
        self.timeMean = {}
        self.timeMed = {}
        self.timeStd = {}

    def getLightCurveData(self, radius=None, return_type='VOTABLE'):
        """
        Class method to get light curve data from a particular source

        Radius of search (optional) set first and then catalogue queried through specific method calls. Will be
        overloaded for different catalogs

        :param radius: Cone search radius in arcseconds. Optional
        :type radius: float
        :param return_type: Type of return format required. Should be 'VOTABLE'
        :type return_type: str
        :return: Successful extract and data ingested
        :rtype: bool
        """
        if radius is None:
            radiusDeg = config.coneRadius
        else:
            radiusDeg = radius / 3600

        if self.archive.name == 'ZTF':
            response = getLightCurveDataZTF(self.AO.pos, radiusDeg)
        elif self.archive.name == 'Pan-STARRS':
            response = getLightCurveDataPanSTARRS(self.AO.pos, radiusDeg, return_type='CSV')
        elif self.archive.name == 'PTF':
            response = getLightCurveDataPTF(self.AO.pos, radiusDeg, return_type)
        else:
            return False

        if response[0]:
            self.table = response[1]
            return True
        else:
            return False

    def set_filtersReturned(self):
        """Determine which filters have returned data for an archive

        Need to use this for """
        self.filtersReturned = self.table[self.archive.filterField].unique().tolist()

    def setID(self):
        for f in self.filtersReturned:
            oidArray = self.table[(self.table[self.archive.filterField] == f)][self.archive.oidField].unique()
            if len(oidArray) != 1:
                maxC = 0
                for id in oidArray:
                    countOID = len(self.table[self.table[self.archive.oidField] == id])
                    if countOID > maxC:
                        self.id[f] = id
                        maxC = countOID
                LClog.info(f'Primary OID {self.id[f]} refmag ({f}mag) used. '
                           f'(Note. Error determining OID due to multiple objects.)')
                # self.table.drop(self.table[(self.table[self.archive.filterField] == f) &
                #                            (self.table[self.archive.oidField] != self.id[f])].index, inplace=True)
            else:
                self.id[f] = oidArray[0]
                LClog.info(f'Primary OID {self.id[f]} refmag ({f}mag) used.')
                # if config.verbose == 'full':
                #     print()

    def setRefMag(self, fileCheck=False):
        # TODO must check in oidlist file if fileCheck is True and save refmag to file

        # oidList = []
        if self.archive.code == 'z':
            for f in self.filtersReturned:
                # oidList.append(str(self.id[f]))

                response = getOIDZTFinfo(str(self.id[f]))
                if len(response):
                    oidData = response
                    # oidRow = oidData[(oidData['oid'] == self.id[f]) & (oidData['filtercode'] == f'z{f}')]
                    self.aRefmag[f] = float(oidData[(oidData['oid'] == self.id[f]) &
                                                    (oidData['filtercode'] == f'z{f}')]['refmag'].values)

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

    def setAltSD(self, col_name, group_col):
        series = self.table.groupby(group_col)[col_name]
        for name, group in series:
            self.altSD[name] = self.median[name] - np.nanpercentile(self.table[
                                                                        (self.table[group_col] == name)
                                                                    ][col_name], 16)

    def setMean(self, col_name, group_col):
        """Method to set the median of data

        Value(s) set within the data structure for each individual filter within the data

        :param group_col:
        :type group_col:
        :param col_name: Column name on which to apply the summary, e.g. 'mag'
        :type col_name: str
        """
        self.mean = self.table.groupby(group_col)[col_name].mean()

    def set3Sigma(self):
        for f in 'gr':
            if self.archive.name == 'ZTF' and f in self.filtersReturned:
                self.sigma3[f] = threeSigma(self.archive, f, self.median[f])

    def analyseTimings(self):
        for f in self.filtersReturned:
            filterData = self.table[self.table[self.archive.filterField] == f].copy()
            filterDataSort = filterData.sort_values(by=[self.archive.timeField])
            timeDiff = filterDataSort[self.archive.timeField].diff()
            if not np.all(np.isnan(timeDiff)):
                self.timeMin[f] = np.nanmin(timeDiff)
                self.timeMax[f] = np.nanmax(timeDiff)
                self.timeMean[f] = np.nanmean(timeDiff)
                self.timeMed[f] = np.nanmedian(timeDiff)
                self.timeStd[f] = np.nanstd(timeDiff)
            else:
                self.timeMin[f] = 0.0
                self.timeMax[f] = 0.0
                self.timeMean[f] = 0.0
                self.timeMed[f] = 0.0
                self.timeStd[f] = 0.0

    def cadenceMin(self, col_name, group_col):
        self.timeMin = self.table.groupby(group_col)[col_name].min()

    def addColourColumn(self, series):
        """Method to add a colour column

        As the data does not have a corresponding colour associated, this is added in a column for
        each row, based on the value of the filter in which the observation was made. This is only used for plots.

        :param series: Column name to use for colour selection
        :type series: str
        """
        c = pd.Series({"g": "green", "r": "red", "i": "indigo", "z": "blue", "y": "black", "R": "orange"})
        self.table['colour'] = self.table[series].map(c)

    def addDateColumn(self, tField):
        self.table['dateTime'] = astroTime.Time(self.table[tField], format='mjd').to_datetime()

    def _getStats(self):
        """Method to set statistics

        Internal method to set statistics for instance of object and archive / catalog. Sets series for filters.

        """
        self.setSamples(self.archive.magField, self.archive.filterField)
        self.setMad(self.archive.magField, self.archive.filterField)
        self.setSD(self.archive.magField, self.archive.filterField)
        self.setMedian(self.archive.magField, self.archive.filterField)
        self.setMean(self.archive.magField, self.archive.filterField)
        self.setAltSD(self.archive.magField, self.archive.filterField)
        self.set3Sigma()
        self.analyseTimings()

    def prepRawData(self):
        self.set_filtersReturned()
        self.setID()
        self.setRefMag()
        self._getStats()

    def getData(self):
        """Method to encapsulate extraction of data

         Data extracted from catalog and summary statistical analysis carried out.
         All data stored in class structure

        :return: Status of data extract
        :rtype: bool
        """

        if self.getLightCurveData():
            self.prepRawData()
            self.dataStatus = True
            return True
        else:
            return False

    def getTable(self):
        return self.table

    def _setFilename(self, mag):
        self.filename = f'data/LC/{mag}/' + \
                        f'{self.archive.name}' + \
                        f'{config.filterSelection}' + \
                        f'{self.AO.objectName}' + \
                        f'.csv'

    def fileExists(self, mag):
        self._setFilename(mag)
        return os.path.isfile(self.filename)

    def saveTable(self, mag):
        if not os.path.exists(f'data/LC/{mag}'):
            os.makedirs(f'data/LC/{mag}')
        self.table.to_csv(self.filename, index=False, header=True)

    def loadTable(self):
        self.table = pd.read_csv(self.filename, header=0)
        self.prepRawData()

    def plot(self, fig, ax, archive, filters):
        """Method to encapsulate the plotting of data

        Sets colour column to distinguish different filters used in data, then sets up plot from given
        X and Y columns. Title is set to object name.

        :param fig: Encapsulating figure object
        :type fig:
        :param ax: Axes object
        :type ax:
        :param archive: Archive object used for configuration of the plot
        :type archive: config.Archive
        """
        if self.dataStatus:  # TODO Need to sort this out for different catalogs
            self.addColourColumn(self.archive.filterField)
            self.addDateColumn(self.archive.timeField)

            for i in filters:
                if len(filters) > 1:
                    filterTable = self.table[self.table[self.archive.filterField] == i]
                    if len(filterTable):
                        ax[filters.index(i)].errorbar(filterTable['dateTime'],
                                                      filterTable[archive.magField],
                                                      filterTable[archive.magErr],
                                                      fmt='none',
                                                      ecolor='grey',
                                                      elinewidth=0.7,
                                                      label=f'{i} error ({archive.name})')
                        ax[filters.index(i)].scatter(filterTable['dateTime'],
                                                     filterTable[archive.magField],
                                                     c=filterTable['colour'],
                                                     marker=archive.marker,
                                                     label=f'{i} filter ({archive.name})')
                        self._plotLines(archive, ax[filters.index(i)], i)

                else:
                    ax.scatter(self.table['dateTime'], self.table[archive.magField],
                               c=self.table['colour'], marker=archive.marker,
                               label=f'{archive.name} {filters} filter')
                    self._plotLines(archive, ax, i)

    def _plotLines(self, archive, ax, fs):
        """Plot horizontal lines

        Insert various horizontal lines on axis passed, based on archive and filter and fill between
        3 sigma if available (currently only ZTF-g)

        :param archive:
        :type archive: config.Archive
        :param ax: axis object
        :type ax: Axes object
        :param fs: filter being plotted
        :type fs: str
        """

        if archive.name == 'ZTF' and fs in 'gr':
            ax.axhline(y=self.median[fs],
                       color='black', linestyle='-', label=f'{archive.name} {fs} filter median')
            self.sigma3 = threeSigma(self.archive, fs, self.median[fs])
            upBound = self.median[fs] + self.sigma3
            lowBound = self.median[fs] - self.sigma3

            ax.axhline(y=upBound, color='grey', linestyle='--', alpha=0.3,
                       label=f'{archive.name} {fs} filter, 3{chr(963)}')
            ax.axhline(y=lowBound, color='grey', linestyle='--', alpha=0.3)
            ax.fill_between(self.table.sort_values(by=['dateTime'])['dateTime'],
                            lowBound, upBound, alpha=0.1)

    def objectOutput(self):
        """Method to encapsulate data output

        Table of summary statistics is sent to console.
        """
        print(f"Archive name: {self.archive.name}")
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
        filterLineOut('Sample interval Min.', self.timeMin)
        filterLineOut('Sample interval Mean', self.timeMean, 2)
        filterLineOut('Sample interval Median', self.timeMed, 2)
        filterLineOut('Sample interval Max.', self.timeMax, 1)
        filterLineOut('Sample interval Std Dev.', self.timeStd, 2)
        print()

    def objectStatSave(self, r):
        """Method to encapsulate summary data save (into table row)
        :param r: AOobject pointer
        :type r:

        """
        r['ZTFObject'] = True
        for f in self.filtersReturned:
            r[f'ZTF{f}oid'] = self.id[f]
            r[f'ZTF{f}filterID'] = f
            r[f'ZTF{f}Samples'] = self.samples[f]
            r[f'ZTF{f}MAD'] = self.mad[f]
            r[f'ZTF{f}SD'] = self.SD[f]
            r[f'ZTF{f}Median'] = self.median[f]
            r[f'ZTF{f}Mean'] = self.mean[f]
            r[f'Stat{f}Included'] = True
            r[f'ZTF{f}RefMag'] = self.aRefmag[f]
            r[f'ZTF{f}med_16pc_diff'] = self.altSD[f]
            r[f'ZTF{f}deltaTMin'] = self.timeMin[f]
            r[f'ZTF{f}deltaTMean'] = self.timeMean[f]
            r[f'ZTF{f}deltaTMedian'] = self.timeMed[f]
            r[f'ZTF{f}deltaTMax'] = self.timeMax[f]
            r[f'ZTF{f}deltaTSD'] = self.timeStd[f]

    def outliersExist(self, a, threshold):
        """Check if archive data has samples outside 3 sigma range

        Function for checking if any sample data returned from archive has detections outside the
        3 sigma range for that mag (as calculated previously from reference data for archive.

        Use threshold values of sample count and mag scaling from config to reject obviously bad data

        :param a: archive
        :type a: config.Archive
        :return: True if samples exist outside the median 3 sigma value
        :rtype: bool
        :param threshold: dict from config
        :type threshold: dict
        """
        response = {}
        self.table['outlier'] = 'inside'
        for f in self.filtersReturned:
            if f == 'i':
                response[f] = False
                continue
            upLim = self.median[f] + self.sigma3[f]
            loLim = self.median[f] - self.sigma3[f]

            self.table.loc[((self.table[a.filterField] == f) &
                            ((self.table[a.magField] - self.table[a.magErr]) >= upLim)),
                           'outlier'] = 'high'

            self.table.loc[((self.table[a.filterField] == f) &
                            ((self.table[a.magField] + self.table[a.magErr]) <= loLim)),
                           'outlier'] = 'low'

            outlierCount = len(self.table[(self.table[a.filterField] == f) & (self.table['outlier'] != 'inside')])

            # check for outliers count between min and max (inclusive) values - initially 0 to 30
            # and also check count is less than a percentage of total samples - initially 10%
            withinCountThreshold = threshold['countMin'] <= outlierCount <= threshold['countMax']
            underTotalCountPC = outlierCount <= int(self.samples[f] * threshold['countPC'])
            if withinCountThreshold and underTotalCountPC:
                response[f] = True
            else:
                if not withinCountThreshold:
                    LClog.info(f'Outlier count ({outlierCount}) for {self.AO.objectName} {f} filter, '
                               f'outside count threshold ({threshold["countMin"]} - {threshold["countMax"]}).')
                else:
                    LClog.info(f'Outlier count ({outlierCount}) for {self.AO.objectName} {f} filter, '
                               f'over count percentage ({threshold["countPC"]}).')
                response[f] = False
        return response

    def johansenTest(self, a):
        pass

    def getMedianRow(self, f):
        a = self.archive
        df = pd.DataFrame(self.table.loc[np.where(self.table['filterID'] == f)])
        df['medDiff'] = abs(df[a.magField] - self.median[f])
        row = df[df['medDiff'] == min(df['medDiff'])].index.item()
        return row

    def getRefImage(self, f):
        oid = self.id[f]
        # get first row for this OID and extract field, etc.
        df = self.table[(self.table['oid'] == oid)]
        field = df.iloc[0]['field']
        filtercode = df.iloc[0]['filtercode']
        ccdid = df.iloc[0]['ccdid']
        quadrant = df.iloc[0]['qid']
        response = getZTFRefImage(field, filtercode, ccdid, quadrant, self.AO.pos, config.imgSize)
        return response

    def getImages(self, row, sizes: list):
        mjd = self.table.loc[row]['mjd']
        filefracday = self.table.loc[row]['filefracday']
        field = self.table.loc[row]['field']
        filtercode = self.table.loc[row]['filtercode']
        ccdid = self.table.loc[row]['ccdid']
        quadrant = self.table.loc[row]['qid']
        response = getZTFImage(mjd, filefracday, field, filtercode, ccdid, quadrant, self.AO.pos, sizes)

        return response
