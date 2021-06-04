"""resample.py

Module
------
resample.py

Description
___________
Module to resample existing SDSS magnitude bins.
Should include

 - load data from all data/LC subdirectories
 - create dataframe
 - sort by g data
 - process stats in sample bins (as per config.stepMag)
 - repeat for r data

Notes
_____
Module to resample existing SDSS magnitude bins to correctly extract and process data from different filters
which will contain data from various magnitudes, dependent upon the 'colour' of the object (star).
e.g. g-r is not a constant, but is indicative of the colour of the object and is used as a proxy to determine
spectral type.
However, since ZTF archive data has been selected based on bins related to the g-mag only but r-mag data also comes
back, this needs to be resampled into separate bins for analysis and output.
This will also allow for changing the mag sample sizes (although obviously won't change the number of samples)

"""
import numpy as np
import pandas as pd
from numpy import arange

from LCExtract import config as LCconfig
from LCExtract.dataretrieve import AstroObjectClass, AODataClass
from SDSSRefs import config
from SDSSRefs.SDSS import SDSSdata, NoLCList, Sky
from SDSSRefs.outputStats import Stats
from SDSSRefs.parse import ParseDir, parseLC, parseSDSS

SDSSfullDataColumns = ['Name', 'RA', 'DEC', 'Description',
                   'SDSSgRefMag', 'SDSSrRefMag',
                   'ZTFObject',
                   'ZTFgoid', 'ZTFgfilterID', 'ZTFgSamples', 'ZTFgMAD', 'ZTFgSD', 'ZTFgMedian', 'ZTFgMean',
                   'StatgIncluded', 'ZTFgRefMag',
                   'ZTFroid', 'ZTFrfilterID', 'ZTFrSamples', 'ZTFrMAD', 'ZTFrSD', 'ZTFrMedian', 'ZTFrMean',
                   'StatrIncluded', 'ZTFrRefMag']

SDSSDataColumns = {'g': ['Name', 'RA', 'DEC', 'Description',
                    'SDSSgRefMag',
                    'ZTFObject',
                    'ZTFgoid', 'ZTFgfilterID', 'ZTFgSamples', 'ZTFgMAD', 'ZTFgSD', 'ZTFgMedian', 'ZTFgMean',
                    'StatgIncluded', 'ZTFgRefMag'],
                   'r': ['Name', 'RA', 'DEC', 'Description',
                    'SDSSrRefMag',
                    'ZTFObject',
                    'ZTFroid', 'ZTFrfilterID', 'ZTFrSamples', 'ZTFrMAD', 'ZTFrSD', 'ZTFrMedian', 'ZTFrMean',
                    'StatrIncluded', 'ZTFrRefMag']}


def getRange(df, magField):
    mod = int(1 / config.stepMag)
    minVal = (int(df[magField].min() * mod)) / mod
    maxVal = (int(df[magField].max() * mod) + 1) / mod
    bins = int((maxVal - minVal) / config.stepMag)
    return minVal, maxVal, bins


def resample():
    statsOfStats = Stats()

    LCconfig.filterSelection = 'gr'
    archive = LCconfig.archives['z']
    skyPlot = Sky()

    if config.LCDataRefresh:
        # Load all LC data into a dataframe
        LCDataTable = pd.DataFrame()
        LCDataTable = parseLC(LCDataTable)
    else:
        SDSSDataTable = pd.DataFrame()
        SDSSDataTable = parseSDSS(SDSSDataTable)  # get SDSS object data

        oidListDF = pd.read_csv(LCconfig.OIDListFile, header=0)  # get ZTF OID data with refmags

    for f in 'gr':
        filterData = SDSSDataTable[(SDSSDataTable[f'ZTF{f}Samples'] != 0)][SDSSDataColumns[f]]  # get data for specific filter
        fullDF = pd.merge(filterData, oidListDF, left_on=[f'ZTF{f}oid'], right_on=['oid'], how='inner')  # merge on OID

        magField = f'ZTF{f}RefMag'
        minMag, maxMag, bins = getRange(fullDF, magField)

        filterData = fullDF[fullDF['filtercode'] == f'z{f}'][magField]
        hist = np.histogram(filterData, bins=bins, range=(minMag, maxMag))

        for mag in arange(minMag, maxMag, config.stepMag):

            lowMagLimit = mag - config.stepMag / 2
            highMagLimit = mag + config.stepMag / 2
            if config.verbose:
                print(f'Sampling {f} filter data, magnitude {mag} from SDSS.')
            filterRange = fullDF[(fullDF[magField].between(lowMagLimit, highMagLimit))]
            if not len(filterRange):
                continue

            sdssDataHolder = SDSSdata(mag)
            sdssDataHolder.setTable(filterRange)
            # status, objectsList = sdssDataHolder.getSamples()
            sdssDataHolder.plotSkyObjects(skyPlot, f)
            # if status:
            #     if config.verbose:
            #         print(f'Sample completed successfully.')
            #     sdssDataHolder.addColumns()
            #     noLClist = NoLCList()
            #     for i in objectsList:
            #         if config.verbose == 'full':
            #             print(f'Querying data for object {i.index + 1} ({i["Name"]}). ', end='')
            #         objectHolder = AstroObjectClass(i['Name'], i['RA'], i['DEC'], i['Description'])
            #         archiveDataHolder = AODataClass(objectHolder, archive)
            #
            #         # retrieve light curve data from file
            #         if config.verbose == 'full':
            #             print(f'\r{" ":66}\r', end='')
            #             print(f'Object {i.index + 1} lightcurve data loaded from previous file.')
            #         archiveDataHolder.loadTable()
            #         archiveDataHolder.objectStatSave(i)

            # sdssDataHolder.saveSamples()
            sdssDataHolder.statsOfStats(f)
            sdssDataHolder.optimiseStats(f)
            statsOfStats.add_row(sdssDataHolder.getStats(f))
            print()

    statsOfStats.removeBlanks()
    statsOfStats.save()
    statsOfStats.plot()
    skyPlot.show()
