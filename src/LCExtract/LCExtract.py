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
from LCExtract.dataretrieve import DataClass
from LCExtract.entry import getObjects, setFilterUsage, setArchiveUsage
from LCExtract.config import archives


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
    for a in archiveList:
        for i in objectsList:
            objectHolder = DataClass(i['Name'], i['RA'], i['DEC'], i['Description'])
            if objectHolder.getData(archives[a]):
                objectHolder.objectOutput(archives[a])
            else:
                print(f'Object name: {i["Name"]} - No data available or retrieved from {archives[a].name}')
                print()
