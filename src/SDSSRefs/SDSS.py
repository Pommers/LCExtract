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
        self.fig.show()


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
        self.totalSamples = 0
        self.medianOfSD = {'g': 0.0, 'r': 0.0}
        self.SDofSD = {'g': 0.0, 'r': 0.0}
        self.meanOfSD = {'g': 0.0, 'r': 0.0}
        self.usedSamples = {'g': 0, 'r': 0}
        self.meanSamplesLC = {'g': 0.0, 'r': 0.0}
        self.ZTFrefMag = {'g': 0.0, 'r': 0.0}
        self.SDSSrefMag = {'g': 0.0, 'r': 0.0}
        self.filename = f'data/SDSS/fullSDSSdata.csv'

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
                f"{self.mag} as Description, g as SDSSgRefMag, r as SDSSgRefMag FROM PhotoPrimary " \
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
            # ascii.write(self.samples,
            #             output=self.filename,
            #             format='csv',
            #             overwrite=True)
        return False

    def statsOfStats(self):
        p = self.samples[(self.samples['ZTFObject'] == True)]

        self.totalSamples = len(p['ZTFObject'])
        for f in 'gr':
            self.usedSamples[f] = len(p[(p[f'ZTF{f}Samples'] != 0)])
            q = p.to_pandas()
            if self.totalSamples and self.usedSamples[f] and not np.all(np.isnan(q[f'ZTF{f}SD'])):
                self.medianOfSD[f] = np.nanmedian(p[f'ZTF{f}SD'])
                self.SDofSD[f] = np.nanstd(p[f'ZTF{f}SD'])
                self.meanOfSD[f] = np.nanmean(p[f'ZTF{f}SD'])
                self.meanSamplesLC[f] = np.nanmean(p[f'ZTF{f}Samples'])

    def optimiseStats(self):
        for f in 'gr':
            if not self.usedSamples[f]:
                continue
            else:
                p = self.samples[(self.samples[f'Stat{f}Included'] == True)]
                # r = p[
                #     (p['ZTFSD'] > (self.meanOfSD - config.sigma * self.SDofSD)) &
                #     (p['ZTFSD'] < (self.meanOfSD + config.sigma * self.SDofSD))
                #     ]
                r = p[p[f'ZTF{f}SD'] < (self.meanOfSD[f] + config.sigma * self.SDofSD[f])]
                if len(r) < self.usedSamples[f]:
                    if len(r):
                        self.usedSamples[f] = len(r)
                        self.medianOfSD[f] = np.nanmedian(p[f'ZTF{f}SD'])
                        self.SDofSD[f] = np.nanstd(r[f'ZTF{f}SD'])
                        self.meanOfSD[f] = np.nanmean(r[f'ZTF{f}SD'])
                        self.meanSamplesLC[f] = np.nanmean(r[f'ZTF{f}Samples'])

    # def optimiseStats_old(self):
    #     while self.totalSamples > config.finalSampleSize:
    #         pass
    #         p = self.samples[(self.samples['StatIncluded'] == True) & (self.samples['ZTFObject'] == True)]
    #         self.samples['StatIncluded'][self.samples['ZTFSD'] == p['ZTFSD'].max()] = False
    #         q = self.samples[(self.samples['StatIncluded'] == True) & (self.samples['ZTFObject'] == True)]
    #
    #         if not len(q):
    #             break
    #         prevSD = self.SDofSD
    #         self.totalSamples = len(q)
    #         self.SDofSD = np.std(q['ZTFSD'])
    #         self.meanSamplesLC = np.mean(q['ZTFSamples'])
    #         if prevSD * config.DQ < self.SDofSD:
    #             if config.verbose:
    #                 print(f'Reached {config.DQ:2.0%} asymptote. ', end='')
    #             break
    #         if q[q['ZTFSD'] < self.SDofSD * 3]:
    #             break
    #     if config.verbose:
    #         print(f'Used {self.totalSamples} samples for magnitude {self.mag:4.2f}±{config.stepMag / 2:5.3f}')
    #     pass

    def getStats(self, f):
        statsList = [self.mag, self.totalSamples]
        for f in 'gr':
            statsList.append(self.medianOfSD[f])
            statsList.append(self.SDofSD[f])
            statsList.append(self.meanOfSD[f])
            statsList.append(self.usedSamples[f])
            statsList.append(self.meanSamplesLC[f])
        return statsList

    def plotSkyObjects(self, sky: Sky):
        ra = coord.Angle(self.samples['RA'] * u.degree)
        ra = ra.wrap_at(180 * u.degree)
        dec = coord.Angle(self.samples['DEC'] * u.degree)  # .filled(np.nan)
        sky.ax.scatter(ra.radian, dec.radian,
                       marker='.',
                       c=self.samples["SDSSgRefMag"],
                       cmap=sky.cmap, alpha=0.9)

    def setTable(self, tableData: pd.DataFrame):
        self.samples = Table.from_pandas(tableData)
