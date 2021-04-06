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
from LCExtract.entry import getObjects, setFilterUsage


def startup():
    print()
    print('Lightcurve data extract\n'
          '-----------------------')
    print()


def LCExtract():

    startup()
    objectsList = getObjects()
    setFilterUsage()
    for i in objectsList:
        objectHolder = DataClass(i['Name'], i['RA'], i['DEC'], i['Description'])
        if objectHolder.getData():
            objectHolder.objectOutput()
        else:
            print(f'Object name: {i["Name"]} - No data available or retrieved')
            print()
