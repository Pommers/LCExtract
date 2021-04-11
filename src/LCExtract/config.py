"""
module
------
config.py

summary
-------
Configuration and global data for package
"""
import collections
from LCExtract.utilities import Namespace

# archAvail = 'pz'  # string for archives available to query

Archive = collections.namedtuple('Archive', 'name code filters URL magField timeField filterField marker')

panstarrs = Archive(name='PanSTARRS', code='p', filters='grizy',
                    URL='https://catalogs.mast.stsci.edu/api/v0.1/panstarrs',
                    magField='psfMag', timeField='obsTime', filterField='filtercode', marker='*')
ztf = Archive(name='ZTF', code='z', filters='gri',
              URL='https://irsa.ipac.caltech.edu/cgi-bin/ZTF/',
              magField='mag', timeField='mjd', filterField='filterID', marker='.')

archives = {'p': panstarrs, 'z': ztf}
archAvail = "".join(list(archives.keys()))

# Global variables
coneRadius = 1 / 3600  # 1 arcseconds
filterSelection = 'grizy'
defaultFileName = 'data/test_new.csv'
badResponse = (False, '')

# baseURL = {'ZTF': 'https://irsa.ipac.caltech.edu/cgi-bin/ZTF/',
#            'PanSTARRS': 'https://catalogs.mast.stsci.edu/api/v0.1/panstarrs'}

# baseURL = Namespace(baseURL)
