"""
module
------
config.py

summary
-------
Configuration and global data for package
"""
import collections

import numpy as np

Archive = collections.namedtuple('Archive', 'name code filters URL magField magErr timeField filterField '
                                            'oidField marker')

# TODO check magnitude error field for PTF and Pan-STARRS (not yet configured)
panstarrs = Archive(name='Pan-STARRS', code='p', filters='grizy',
                    URL='https://catalogs.mast.stsci.edu/api/v0.1/panstarrs',
                    magField='psfMag', magErr='magerr', timeField='obsTime',
                    filterField='filtercode', oidField='objID', marker='*')
ztf = Archive(name='ZTF', code='z', filters='gri',
              URL='https://irsa.ipac.caltech.edu/cgi-bin/ZTF/',
              magField='mag', magErr='magerr', timeField='mjd',
              filterField='filterID', oidField='oid', marker='.')
ptf = Archive(name='PTF', code='r', filters='gR',
              URL='https://irsa.ipac.caltech.edu/cgi-bin/Gator/',
              magField='mag_autocorr', magErr='magerr', timeField='obsmjd',
              filterField='filterID', oidField='oid', marker='+')

archives = {'p': panstarrs, 'z': ztf, 'r': ptf}
archAvail = "".join(list(archives.keys()))

# Global variables
coneRadius = 1 / 3600  # 1 arcseconds
defaultFileName = 'data/test_objects.csv'
badResponse = (False, '')

# filter management
filterSelection = 'grizyR'  # filters selected for run
filtersAvail = {}  # filters determined as subset of selection and archive
filtersRet = {}  # filters actually returned per archive (may not be used)

LCDataSave = True
LCDataRefresh = True

OIDListFile = 'data/LC/oidList.csv'

# median SD of SDSS point source (LSfit 8th order poly coefficients)
medSD_c = (2.19790659e-07, -2.89677718e-05,  1.65032595e-03, -5.30759586e-02,
           1.05382668e+00, -1.32260634e+01,  1.02451798e+02, -4.47763281e+02, 8.45226949e+02)

# set threshold values for determining detections of interest in lightcurves
threshold = {'countMin': 1, 'countMax': 30, 'countPC': 0.1, 'magScalar': 1.0}

# SDSS configuration
sampleSize = 20
magWindow = 0.5  # ± mag
timeWindow = 0.5  # ± days
regionSize = 5 / 60  # search box ± in arcmin
