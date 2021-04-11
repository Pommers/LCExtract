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

# Global variables
archAvail = 'pz'  # string for archives available to query

Archive = collections.namedtuple('Archive', 'name code filters URL magField timeField')

panstarrs = Archive(name='PanSTARRS', code='p', filters='grizy',
                    URL='https://catalogs.mast.stsci.edu/api/v0.1/panstarrs',
                    magField='mag', timeField='obsTime')
ztf = Archive(name='ZTF', code='z', filters='gri',
              URL='https://irsa.ipac.caltech.edu/cgi-bin/ZTF/',
              magField='mag', timeField='mjd')

archives = {'p': panstarrs, 'z': ztf}

filterSelection = 'grizy'
defaultFileName = 'data/test_new.csv'

baseURL = {'ZTF': 'https://irsa.ipac.caltech.edu/cgi-bin/ZTF/',
           'PanSTARRS': 'https://catalogs.mast.stsci.edu/api/v0.1/panstarrs'}

baseURL = Namespace(baseURL)
