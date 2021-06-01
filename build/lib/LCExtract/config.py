"""
module
------
config.py

summary
-------
Configuration and global data for package
"""
import collections

# archAvail = 'pz'  # string for archives available to query

Archive = collections.namedtuple('Archive', 'name code filters URL magField timeField filterField marker')

panstarrs = Archive(name='Pan-STARRS', code='p', filters='grizy',
                    URL='https://catalogs.mast.stsci.edu/api/v0.1/panstarrs',
                    magField='psfMag', timeField='obsTime', filterField='filtercode', marker='*')
ztf = Archive(name='ZTF', code='z', filters='gri',
              URL='https://irsa.ipac.caltech.edu/cgi-bin/ZTF/',
              magField='mag', timeField='mjd', filterField='filterID', marker='.')
ptf = Archive(name='PTF', code='r', filters='gR',
              URL='https://irsa.ipac.caltech.edu/cgi-bin/Gator/',
              magField='mag_autocorr', timeField='obsmjd', filterField='filterID', marker='+')

archives = {'p': panstarrs, 'z': ztf, 'r': ptf}
archAvail = "".join(list(archives.keys()))

# Global variables
coneRadius = 1 / 3600  # 1 arcseconds
filterSelection = 'grizyR'
defaultFileName = 'data/test_new.csv'
badResponse = (False, '')

# baseURL = {'ZTF': 'https://irsa.ipac.caltech.edu/cgi-bin/ZTF/',
#            'PanSTARRS': 'https://catalogs.mast.stsci.edu/api/v0.1/panstarrs'}

# baseURL = Namespace(baseURL)
