"""
Module
------
entry.py:

Summary
-------
Contains routines to facilitate entry of required input data objects, either from manual user input or from
a data file (csv).

Notes
-----
May include database access in the future to save having to create a csv from catalogue after update

"""
# Imports included here:
import re

# Astropy
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates.name_resolve import NameResolveError
from astropy.io import ascii
from astropy.table import Table

from LCExtract import config


def setArchiveUsage():
    print('Please specify which archives to query.')
    for a in config.archives:
        print(f'{config.archives[a].code} - {config.archives[a].name}')
    getch = input(f'Enter for default ({config.archAvail})')
    archives = config.archAvail if getch == '' else re.findall(f'[{config.archAvail}]', getch.lower())
    print()
    return archives


def setFilterUsage():
    # uses global filterSelection

    getch = input(f'Please select filters to display (deafult is {config.filterSelection})....: ')
    # if no entry set use default for session
    if getch != '':
        config.filterSelection = getch
    print()


def setEntryType():
    print('Script will accept file or manual input. Default is manual.')
    while True:
        getch = input(f'Please select file (f) or manual object (m) entry..........: ')
        if getch == '':
            getch = 'm'
        if getch[0].lower() in ('f', 'm'):
            break
    print()

    return getch[0].lower()


def setManualEntryType():
    while True:
        getch = input(f'Please select named object (n) or coordinate (c) entry..........: ')
        if getch == '':
            continue
        if getch[0].lower() in ('n', 'c'):
            break
    print()

    return getch[0].lower()


def getObjectsCSV():
    # uses global defaultFileName
    error_to_catch = getattr(__builtins__, 'FileNotFoundError', IOError)
    while True:
        getch = input(f'Please enter filename, or <CR> for default ({config.defaultFileName})...: ')
        getch = config.defaultFileName if getch == '' else getch
        try:
            f = open(getch)
        except error_to_catch:
            print(f'Unable to locate file "{getch}". Please try again.')
        else:
            f.close()
            print(f'Using file "{getch}".')
            print()
            break

    data = ascii.read(getch, guess=False, format='csv', header_start=0, data_start=1)
    return data


def getUserObject():
    manualEntryType = setManualEntryType()
    if manualEntryType == 'n':  # named object entry
        while True:
            tempName = input('Enter object name....: ')
            try:
                c = SkyCoord.from_name(tempName, parse=True)
            except NameResolveError:
                print('Unable to resolve. Please try again.')
            else:
                print(f'Object {tempName} found in catalog.')
                break
    elif manualEntryType == 'c':  # object coordinate entry
        while True:
            print('Enter object coordinates (ICRS frame. Deg assumed unless specified).')
            tempCoordRA = input('RA (e.g. 10.625, 10d37m30s, 0h42m30s, 00 42 30)....: ')
            tempCoordRA += 'd' if not re.findall('[hdms]', tempCoordRA) else ''
            tempCoordDEC = input('DEC (e.g. 41.2, 41d12m00s, +41 12 00)...: ')
            tempCoordDEC += 'd' if not re.findall('[dms]', tempCoordDEC) else ''
            try:
                c = SkyCoord(tempCoordRA, tempCoordDEC)
            except ValueError:
                print('Unable to identify position. Please try again.')
            except u.UnitsError:
                print('Units error occurred. Please try again.')
            else:
                print(f'Object at {c.to_string("hmsdms")} found.')
                tempName = 'Object in ' + c.get_constellation()
                break
    print()
    manual = [{'Name': tempName,  # 'Name': 'Sky position: 153.139, 53.117',
               'RA': c.ra.degree,  # 'RA': 153.1393271,
               'DEC': c.dec.degree,  # 'DEC': 53.117343,
               'Description': f'Position: {c.to_string("hmsdms")} '}]  # 'Description': 'Manual input test'
    return manual


def getObjects():
    entryType = setEntryType()
    if entryType == 'f':
        return getObjectsCSV()
    elif entryType == 'm':
        singleObjectData = getUserObject()

        tbl = Table(rows=singleObjectData)
        return tbl
