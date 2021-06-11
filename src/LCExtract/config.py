"""
module
------
config.py

summary
-------
Configuration and global data for package
"""
import collections
import logging
import os

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
defaultFileName = 'data/test_single.csv'
badResponse = (False, '')

# filter management
filterSelection = 'grizyR'  # filters selected for run
filtersAvail = {}  # filters determined as subset of selection and archive
filtersRet = {}  # filters actually returned per archive (may not be used)

LCDataSave = True
LCDataRefresh = True

OIDListFile = 'data/LC/oidList.csv'

# median SD of SDSS point source (LSfit 8th order poly coefficients)
# medSD_c = [2.19790659e-07, -2.89677718e-05, 1.65032595e-03, -5.30759586e-02,
#            1.05382668e+00, -1.32260634e+01, 1.02451798e+02, -4.47763281e+02, 8.45226949e+02]

medSD_c = {'g': (5.52593602e-07, -7.84244669e-05,  4.84231576e-03, -1.69929569e-01,
                 3.70746577e+00, -5.15010235e+01,  4.44837098e+02, -2.18429205e+03,
                 4.66813446e+03),
           'r': (-5.42286684e-08,  7.77543805e-06, -4.84682328e-04,  1.70914304e-02,
                 -3.71645230e-01,  5.08729087e+00, -4.26875655e+01,  2.00139505e+02,
                 -3.99936057e+02)
           }

# set threshold values for determining detections of interest in lightcurves
threshold = {'countMin': 1, 'countMax': 30, 'countPC': 0.1, 'magScalar': 1.0}

# SDSS configuration
sampleSize = 20
magWindow = 0.5  # ± mag
timeWindow = 0.5  # ± days
regionSize = 5 / 60  # search box ± in arcmin

#  verbose setting for output monitoring and debug
verbose = 'full'  # use 'full' / 'minimal' / False

# setup logging for LCExtract
if not os.path.exists('logs'):
    os.makedirs('logs')
LClog = logging.getLogger(__name__)
LClog.setLevel(logging.DEBUG)  # DEBUG / INFO / WARNING / ERROR / CRITICAL

formatter = logging.Formatter('%(asctime)s: %(levelname)s: %(module)s: %(message)s')
file_handler = logging.FileHandler('logs/LCExtract.log')
file_handler.setLevel(logging.INFO)
file_handler.setFormatter(formatter)

stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)

LClog.addHandler(file_handler)
LClog.addHandler(stream_handler)
