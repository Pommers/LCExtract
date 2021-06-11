"""
Module
------
LCExtract.py: Wrapper for lightcurve data extraction package

Summary
-------
    Allows selection of file or manual input.

    Allows filter selection.

    Iterates through required entries

    retrieving available data

    summarising and plotting output
"""
# Self-authored package modules for inclusion
from LCExtract import config
from LCExtract.config import LClog
from LCExtract.SDSS import SDSSdata
from LCExtract.config import archives
from LCExtract.dataretrieve import AstroObjectClass, AODataClass
from LCExtract.entry import getObjects, setFilterUsage, setArchiveUsage
from LCExtract.ztf import refZTFobj, getOIDZTFinfo, getZTFOidOffsetStats


def startup():
    print()
    print('Lightcurve data extract\n'
          '-----------------------')
    print()


def LCExtract():
    archiveLCData = {}
    outliersFound = {}
    refAOList = {}
    meanMagOffset = {}
    startup()
    objectsList = getObjects()
    archiveList = setArchiveUsage()
    setFilterUsage()
    for i in objectsList:
        AO = AstroObjectClass(i['Name'], i['RA'], i['DEC'], i['Description'])
        print(f"\n-------------------------------------------------------------\n"
              f"Object name: {AO.objectName} - summary statistics")

        for a in archiveList:
            archiveLCData[a] = AODataClass(AO, archives[a])
            if archiveLCData[a].getData():
                AO.incrPlotRows()
                AO.append_filtersToPlot(archiveLCData[a].filtersReturned)
                archiveLCData[a].objectOutput()
                if archives[a].name == 'ZTF':
                    # check if archive data has data outside 3 sigma
                    outliersFound[a] = archiveLCData[a].outliersExist(archives[a], config.threshold)
                    if outliersFound[a]:
                        outliers = archiveLCData[a].table[archiveLCData[a].table['outlier'] != 'inside']
                        # for each of the outlier datapoints in the LC sample
                        for i1, s in outliers.iterrows():
                            refs = refZTFobj(AO.pos,
                                             s[archives[a].timeField],
                                             s['filefracday'],
                                             s[archives[a].magField])
                            ZTFoffset, ZTFsd = getZTFOidOffsetStats(refs)
                            print(f'{AO.objectName} - Outlier MJD{s[archives[a].timeField]}: '
                                  f'Offset = {ZTFoffset:-8.5f} '
                                  f'SD = {ZTFsd:-8.5f}')
                        sdssRefAOs = SDSSdata(archiveLCData[a].median, AO.pos)
                        if sdssRefAOs.getSamples():
                            pass
                        # return a df of image data frames (idf) to check
                        # for idf in archiveLCData[a].outlierIDFList:
                        # refAOList[idf] = archiveLCData[a].getRefObjects(idf)
                        # get a df of ref objects in idf (with SDSS mag)
                        # meanMagOffset[idf] = archiveLCData[a].getIDFmagOffset(refAOList[idf])
                        # get mean mag offset for ref objects in IDF
            else:
                print(f'No data available or retrieved from {archives[a].name}')
                print()
        filtersToPlot = AO.get_filtersInAO()
        if AO.plotRowsTrue():
            fig, ax = AO.preparePlot(filtersToPlot)
            for a in archiveList:
                archiveLCData[a].plot(fig, ax, archives[a], filtersToPlot)
            AO.finalisePlot(fig, ax, filtersToPlot)
