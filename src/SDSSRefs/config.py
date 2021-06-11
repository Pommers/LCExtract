"""config.py

Module
------
config.py

Description
-----------
Module for SDSS reference sample extraction to contain configuration data values

"""
import logging
from datetime import date
import os

# working directory for data
wd = f"data{os.sep}"

# main method to execute
parse = False

# no. of samples
startingSampleSize = 200
finalSampleSize = 10

# SDSS data storage
SDSSDataSave = True
SDSSDataRefresh = True
SDSSBasePath = f'{wd}SDSS'
SDSSDataFilename = f'{SDSSBasePath}{os.sep}fullSDSSdata_{startingSampleSize:03}.csv'

# ZTF Scatter plot filename
plotsPath = f'{wd}plots'
ZTFScFilename = f'{plotsPath}{os.sep}ztf_scatter_{startingSampleSize:03}_gr.png'
ZTFSc2Filename = f'{plotsPath}{os.sep}ztf_scatter2_{startingSampleSize:03}_gr.png'


# Sky plot filename
skyplotFilename = f'{plotsPath}{os.sep}skyplot_{startingSampleSize:03}_gr.png'

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
optCycles = 3

# data quality
sigma = 3  # normal distribution sample check fit

# LC data storage
LCDataSave = True
LCDataRefresh = False
LCNotFoundFile = 'data/LC/NoLCList.csv'

#  verbose setting for output monitoring and debug
verbose = 'minimal'  # use 'full' / 'minimal' / False

# setup logging for SDSSRefs
if not os.path.exists('logs'):
    os.makedirs('logs')
SDSSlog = logging.getLogger(__name__)
SDSSlog.setLevel(logging.INFO)  # DEBUG / INFO / WARNING / ERROR / CRITICAL

formatter = logging.Formatter('%(asctime)s: %(levelname)s: %(module)s: %(message)s')
file_handler = logging.FileHandler('logs/SDSSRefs.log')
file_handler.setLevel(logging.INFO)
file_handler.setFormatter(formatter)

stream_handler = logging.StreamHandler()
stream_handler.setLevel(logging.DEBUG)
stream_handler.setFormatter(formatter)

SDSSlog.addHandler(file_handler)
SDSSlog.addHandler(stream_handler)


# file structure for stats save - suggest date.archive.maxsamples.stepmag
today = date.today()
filename = f"{wd}{today.strftime('%Y%m%d')}.ZTF.{startingSampleSize}.{stepMag}.{sigma}s.csv"
