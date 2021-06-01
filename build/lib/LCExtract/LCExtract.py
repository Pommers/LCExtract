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
from LCExtract.config import archives
from LCExtract.dataretrieve import AstroObjectClass, AODataClass
from LCExtract.entry import getObjects, setFilterUsage, setArchiveUsage


def startup():
    print()
    print('Lightcurve data extract\n'
          '-----------------------')
    print()


def LCExtract():
    startup()
    objectsList = getObjects()
    archiveList = setArchiveUsage()
    setFilterUsage()
    for i in objectsList:
        objectHolder = AstroObjectClass(i['Name'], i['RA'], i['DEC'], i['Description'])
        print(f"Object name: {objectHolder.objectName} - summary statistics")

        plotData = False
        fig, ax = objectHolder.preparePlot(len(config.filterSelection))
        for a in archiveList:
            dataHolder = AODataClass(objectHolder)
            if dataHolder.getData(archives[a]):
                plotData = True
                dataHolder.objectOutput(archives[a])
                dataHolder.getTable()
                dataHolder.plot(fig, ax, archives[a])
            else:
                print(f'No data available or retrieved from {archives[a].name}')
                print()

        objectHolder.finalisePlot(plotData, fig, ax)
