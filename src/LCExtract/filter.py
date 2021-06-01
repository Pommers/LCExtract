"""filter.py

Module
------
filter.py

Description
-----------
Module for the management of filter requests within the lightcurve extraction and output.
Each archive / catalogue will have its own filter availability, and the user may select all
or a subset of these. In addition, not all lightcurve extracts / returns will contain samples
with all filters included. Plotting and output of data in different filters is thus required.
"""
import re
import numpy as np

from LCExtract import config


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
