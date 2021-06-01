"""config.py

Module
------
config.py

Description
-----------
Module for SDSS reference sample extraction to contain configuration data values

"""
from datetime import date


# main method to execute
parse = False

# SDSS data storage
SDSSDataSave = True
SDSSDataRefresh = True

# no. of samples
startingSampleSize = 200
finalSampleSize = 10

# RA limits (probably won't be used)
RAmin = 0.0
RAmax = 360.0

# Dec limits (conservative limits for observing for ZTF)
DecMin = -20.0
DecMax = 70.0

# min/max magnitude
minMag = 12
maxMag = 25.0

# step size
stepMag = 0.25

# iterations (to get to a min)
cycleIter = 4

# data quality
sigma = 3  # normal distribution sample check fit

# LC data storage
LCDataSave = True
LCDataRefresh = False
LCNotFoundFile = 'data/LC/NoLCList.csv'

#  verbose setting for output monitoring and debug
verbose = 'full'  # use 'full' / 'minimal' / False

# file structure for stats save - suggest date.archive.maxsamples.stepmag
today = date.today()
wd = "data/"
filename = f"{wd}{today.strftime('%Y%m%d')}.ZTF.{startingSampleSize}.{stepMag}.{sigma}s.csv"
