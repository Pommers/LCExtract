"""SDSS.py

Module
------
SDSS.py

Description
-----------
Module for access to SDSS data. Contains SDSSdata Class

"""
import math
import os

import numpy as np
import pandas as pd
from astropy import coordinates as coord
from astropy import units as u
from astropy.io import ascii
from astropy.table import Table
from astroquery.sdss import SDSS
from matplotlib import cm
from matplotlib import colors
from matplotlib import pyplot as plt

from SDSSRefs import config
from LCExtract import config as LCconfig


class Sky:
    def __init__(self):
        self.cmap = cm.get_cmap('gist_yarg')
        self.fig, self.ax = self.create()

    def create(self):
        self.fig = plt.figure(figsize=(10, 5), dpi=300)
        self.ax = self.fig.add_subplot(111, projection="mollweide")
        self.ax.set_title("SDSS Object location")
        self.ax.set_xticklabels(['14h', '16h', '18h', '20h', '22h', '0h', '2h', '4h', '6h', '8h', '10h'])
        self.ax.grid(True)
        norm = colors.Normalize(vmin=11.5, vmax=24.5)
        self.fig.colorbar(cm.ScalarMappable(norm=norm, cmap=self.cmap), label='g/r [mag]', shrink=.6, pad=.05,
                          aspect=15)
        return self.fig, self.ax

    def ax(self):
        return self.ax

    def show(self):
        if not os.path.isdir(config.plotsPath):
            os.makedirs(config.plotsPath)

        self.fig.show()
        self.fig.savefig(config.skyplotFilename)


class Plot2:
    def __init__(self):
        self.fig, self.ax = self.create()

    def create(self):
        if not os.path.isdir(config.plotsPath):
            os.makedirs(config.plotsPath)

        self.fig, self.ax = plt.subplots(2, 1, sharex=True, figsize=(8, 7), dpi=300)

        return self.fig, self.ax

    def filterPlot(self, f, filterDF, stats):
        i = 'gr'.index(f)
        filterData = stats[stats['filter'] == f]
        # plot Scatter of ZTF lightcurve samples based on SD of SD

        # mc = np.where((filterDF[f'Stat{f}Included']), f, 'black')
        inc = filterDF[(filterDF[f'Stat{f}Included'] == True)]
        exc = filterDF[(filterDF[f'Stat{f}Included'] == False)]
        # filterData[filterData['ZTFgfilterID'] == 'g']['ZTFgRefMag']['ZTFgSD']['ZTFgSamples']
        self.ax[i].scatter(inc[f'ZTF{f}RefMag'], inc[f'ZTF{f}SD'], label='Light curve object SDs',
                           alpha=0.4, marker='.', c=f)
        self.ax[i].scatter(exc[f'ZTF{f}RefMag'], exc[f'ZTF{f}SD'],
                           label=f'Excluded Light curves (>3{chr(963)})',
                           alpha=0.4, marker='.', c='black')

        upBound = filterData[f'MedianOfSD'] + filterData[f'SDofSD']
        lowBound = filterData[f'MedianOfSD'] - filterData[f'SDofSD']
        self.ax[i].fill_between(filterData['mag'], upBound, lowBound, color='grey', alpha=0.7,
                                label=f'{f}-band sample {chr(963)}')
        self.ax[i].plot(filterData['mag'], filterData[f'MedianOfSD'],
                        label=f'Median of {f}-band', c='black', marker='.')
        self.ax[i].set_ylabel(f'Std Dev ({chr(963)})', fontsize=10)
        self.ax[i].set_title('ZTF Lightcurves Std Dev', fontsize=12)
        self.ax[i].legend(loc='upper left')
        self.ax[i].set_ylim(bottom=0, top=max(upBound) * 1.1)
        # bestFit(filterData['mag'], filterData[f'MedianOfSD'], self.ax[i])

    def show(self):
        self.fig.tight_layout()
        self.fig.show()
        self.fig.savefig(config.ZTFSc2Filename)


class NoLCList:
    def __init__(self):
        self.file = self.openNLC()
        self.list = Table
        self.getFile()

    @staticmethod
    def openNLC():
        try:
            f1 = open(config.LCNotFoundFile, 'a')
        except FileNotFoundError:
            if not os.path.exists('data/LC'):
                os.makedirs('data/LC')
            if not os.path.isfile(config.LCNotFoundFile):
                f1 = open(config.LCNotFoundFile, 'a')
                f1.write('ID,mag' + os.linesep)
                f1.flush()
        return f1

    def getFile(self):
        self.list = ascii.read(config.LCNotFoundFile, format='csv')

    def add(self, ID, mag):
        self.file.write(f'{ID},{mag}' + os.linesep)

    def inList(self, ID):
        find = (self.list['ID'] == ID)
        return len(self.list[find]) > 0


class SDSSdata:
    def __init__(self, mag):

        self.samples = Table
        self.mag = mag
        self.filter = ''
        self.totalSamples = {'g': 0, 'r': 0}
        self.medianOfSD = {'g': 0.0, 'r': 0.0}
        self.SDofSD = {'g': 0.0, 'r': 0.0}
        self.meanOfSD = {'g': 0.0, 'r': 0.0}
        self.usedSamples = {'g': 0, 'r': 0}
        self.meanSamplesLC = {'g': 0.0, 'r': 0.0}
        self.ZTFrefMag = {'g': 0.0, 'r': 0.0}
        self.SDSSrefMag = {'g': 0.0, 'r': 0.0}
        self.filename = config.SDSSDataFilename

    def addColumns(self):
        self.samples.add_column(False, name='ZTFObject')
        for f in 'gr':
            self.samples.add_column(0, name=f'ZTF{f}oid')
            self.samples.add_column('', name=f'ZTF{f}filterID')
            self.samples.add_column(0, name=f'ZTF{f}Samples')
            self.samples.add_column(0.0, name=f'ZTF{f}MAD')
            self.samples.add_column(0.0, name=f'ZTF{f}SD')
            self.samples.add_column(0.0, name=f'ZTF{f}Median')
            self.samples.add_column(0.0, name=f'ZTF{f}Mean')
            self.samples.add_column(False, name=f'Stat{f}Included')
            self.samples.add_column(0.0, name=f'ZTF{f}RefMag')

    def getSamples(self):
        sampleCount = config.startingSampleSize
        stepSize = config.stepMag
        magLow = self.mag - config.stepMag / 2
        magHigh = self.mag + config.stepMag / 2
        query = f"SELECT TOP {sampleCount} objid as Name, ra as RA, dec as DEC, " \
                f"{self.mag} as Description, g as SDSSgRefMag, r as SDSSrRefMag FROM PhotoPrimary " \
                f"WHERE g BETWEEN {magLow} AND {magHigh} " \
                f"AND dec BETWEEN {config.DecMin} AND {config.DecMax} " \
                f"AND type = 6 " \
                f"AND mode = 1"
        self.samples = SDSS.query_sql(query, data_release=13)
        # self.saveSamples()
        return len(self.samples) == sampleCount, self.samples

    def createFile(self, overwrite=False):
        if not os.path.exists('data/SDSS'):
            os.makedirs('data/SDSS')
        # self.fileHandle = open(self.filename, mode='a')

    def saveSamples(self, initialise):
        if initialise:
            self.samples.write(self.filename, format='ascii')  # write first data with header (I hope)
        else:
            with open(self.filename, mode='a') as f:
                f.seek(0, os.SEEK_END)
                self.samples.write(f, format='ascii.no_header')
        return False

    def _setDetailStats(self, p, f):
        self.medianOfSD[f] = np.nanmedian(p[f'ZTF{f}SD'])
        self.SDofSD[f] = np.nanstd(p[f'ZTF{f}SD'])
        self.meanOfSD[f] = np.nanmean(p[f'ZTF{f}SD'])
        self.meanSamplesLC[f] = np.nanmean(p[f'ZTF{f}Samples'])

    def statsOfStats(self, f):
        self.samples[f'Stat{f}Included'][np.where(np.isnan(self.samples[f'ZTF{f}SD']))] = False
        p = self.samples[(self.samples['ZTFObject'] == True)]

        self.totalSamples[f] = len(p['ZTFObject'])

        self.usedSamples[f] = len(p[(p[f'ZTF{f}Samples'] != 0)])
        q = p.to_pandas()
        if self.totalSamples[f] and self.usedSamples[f] and not np.all(np.isnan(q[f'ZTF{f}SD'])):
            self._setDetailStats(p, f)

    def optimiseStats(self, f):
        if self.usedSamples[f]:
            p = self.samples[(self.samples[f'Stat{f}Included'] == True)]
            count = 0
            while len(p) and count < config.optCycles:
                self.samples[f'Stat{f}Included'][
                    np.where(
                        (self.samples[f'ZTF{f}SD'] > (self.medianOfSD[f] + config.sigma * self.SDofSD[f])) |
                        (self.samples[f'ZTF{f}SD'] < (self.medianOfSD[f] - config.sigma * self.SDofSD[f])))
                ] = False
                p = self.samples[(self.samples[f'Stat{f}Included'] == True)]

                # r = p[p[f'ZTF{f}SD'] < (self.meanOfSD[f] + config.sigma * self.SDofSD[f])]
                if len(p) and len(p) != self.usedSamples[f]:
                    self.usedSamples[f] = len(p)
                    self._setDetailStats(p, f)
                else:
                    break
                count += 1

    def getStats(self, f):
        statsList = [self.mag, f, self.totalSamples[f], self.medianOfSD[f], self.SDofSD[f], self.meanOfSD[f],
                     self.usedSamples[f], self.meanSamplesLC[f]]
        return statsList

    def plotSkyObjects(self, sky: Sky, f):
        ra = coord.Angle(self.samples['RA'] * u.degree)
        ra = ra.wrap_at(180 * u.degree)
        dec = coord.Angle(self.samples['DEC'] * u.degree)  # .filled(np.nan)
        sky.ax.scatter(ra.radian, dec.radian,
                       marker='.',
                       c=self.samples[f"SDSS{f}RefMag"],
                       cmap=sky.cmap, alpha=0.9)

    def setTable(self, tableData: pd.DataFrame):
        self.samples = Table.from_pandas(tableData)
