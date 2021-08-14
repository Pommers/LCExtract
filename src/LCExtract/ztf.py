import io
import os
from urllib.error import HTTPError
from urllib.request import urlopen

import astropy.time
import pandas as pd
import pyvo
import requests
from astropy.io import fits

from LCExtract import config
from LCExtract.coord import CoordClass, to_string
from LCExtract.filter import getFilterStr
from LCExtract.utilities import Spinner

delim = "%20"


class OIDListFile:
    header = ['cntr', 'oid', 'ra', 'dec', 'htm20', 'field', 'ccdid', 'qid', 'fid', 'filtercode',
              'x', 'y', 'z', 'ngoodobs', 'ngoodobsrel', 'nobs', 'nobsrel', 'refchi', 'refmag',
              'refmagerr', 'refsharp', 'refsnr', 'astrometricrms', 'chisq', 'con', 'lineartrend',
              'magrms', 'maxmag', 'maxslope', 'meanmag', 'medianabsdev', 'medianmag', 'medmagerr',
              'minmag', 'nabovemeanbystd_1', 'nabovemeanbystd_3', 'nabovemeanbystd_5',
              'nbelowmeanbystd_1', 'nbelowmeanbystd_3', 'nbelowmeanbystd_5', 'nconsecabovemeanbystd_1',
              'nconsecabovemeanbystd_3', 'nconsecabovemeanbystd_5', 'nconsecbelowmeanbystd_1',
              'nconsecbelowmeanbystd_3', 'nconsecbelowmeanbystd_5', 'nconsecfrommeanbystd_1',
              'nconsecfrommeanbystd_3', 'nconsecfrommeanbystd_5', 'nmedianbufferrange',
              'npairposslope', 'percentiles_05', 'percentiles_10', 'percentiles_175',
              'percentiles_25', 'percentiles_325', 'percentiles_40', 'percentiles_60',
              'percentiles_675', 'percentiles_75', 'percentiles_825', 'percentiles_90',
              'percentiles_95', 'skewness', 'smallkurtosis', 'stetsonj', 'stetsonk',
              'vonneumannratio', 'weightedmagrms', 'weightedmeanmag']

    def __init__(self):
        self.openOIDlist()
        self.df = pd.DataFrame()
        self.getFile()

    @staticmethod
    def openOIDlist():
        try:
            f2 = open(config.OIDListFile, 'r')
        except FileNotFoundError:
            if not os.path.exists('data/LC'):
                os.makedirs('data/LC')
            if not os.path.isfile(config.OIDListFile):
                with open(config.OIDListFile, 'a') as f:
                    pd.DataFrame(columns=OIDListFile.header,
                                 index=None).to_csv(f,
                                                    index=False,
                                                    header=OIDListFile.header)  # header=f.tell() == 0

    def getFile(self):
        self.df = pd.read_csv(config.OIDListFile, header=0)

    def add(self, response):
        with open(config.OIDListFile, 'a') as f:
            pd.DataFrame(response).to_csv(f, index=False, header=False,
                                          columns=OIDListFile.header)  # header=f.tell() == 0

    def inList(self, OID):
        if len(self.df):
            find = (self.df['oid'] == int(OID))
        else:
            find = self.df[self.df['oid'] == int(OID)]
        return self.df[find]


oidListFile = OIDListFile()


def getLightCurveDataZTF(coordinates: CoordClass, radius, column_filters=None):
    """Zwicky Transient facility light curve data retrieval

    IRSA provides access to the ZTF collection of lightcurve data through an application program interface (API).
    Search, restriction, and formatting parameters are specified in an HTTP URL. The output is a table in the
    requested format containing lightcurve data satisfying the search constraints.

    Ref. https://irsa.ipac.caltech.edu/docs/program_interface/ztf_lightcurve_api.html

    :param coordinates: Coordinates of object expressed CoordClass notation in J2000 RA Dec (Decimal) format.
    :type coordinates: CoordClass
    :param radius: Radius of cone search ** in degrees ** for passing to ZTF
    :type radius: float
    :param column_filters: Not used currently
    :returns:
        (boolean) Valid data return
        (DataFrame) Data payload
    :rtype: tuple

    """
    filterStr = getFilterStr(config.ztf.filters)  # limit filters (requested) to ZTF subset

    if not filterStr:  # Check table actually has data in it (i.e. possible no lightcurve data exists)
        return config.badResponse

    status = True

    ra = coordinates.ra_str() + delim
    dec = coordinates.dec_str() + delim
    radius_str = to_string(radius, 5)
    # if column_filters is None:
    #     column_filters = {}

    queryPart = "nph_light_curves"
    pos = "POS=CIRCLE" + delim + ra + dec + radius_str
    bandName = "BANDNAME=" + filterStr
    form = "FORMAT=CSV"  # was VOTABLE but changed to CSV
    badCatFlagsMask = "BAD_CATFLAGS_MASK=32768"

    url_payload = f"{config.ztf.URL}{queryPart}?{pos}&{bandName}&{form}&{badCatFlagsMask}"

    if column_filters is not None:
        url_payload += '&' + '&'.join(column_filters)

    # establish http connection
    # http = urllib3.PoolManager()
    # siteData = http.request('GET', url_payload)
    print('Requesting data from Zwicky Transient Facility. Please wait ... ', end='')
    with Spinner():
        try:
            siteData = urlopen(url_payload)
            print(f'\r{" ":66}\r ', end='')
            # print(' ', end='')
        except HTTPError as err:
            if err.code == 400:
                print('Sorry. Could not complete request.')
                return config.badResponse
            else:
                return config.badResponse

    if siteData.status != 200:  # Ensure good response is received back from IRSA
        return config.badResponse

    # read CSV data returned from URL, decode from byte to string format, place in a wrapper and read into DF
    tablePD = pd.read_csv(io.StringIO(siteData.read().decode()))

    if not len(tablePD):  # Check table actually has data in it (i.e. possible no lightcurve data exists)
        return config.badResponse

    # can't remove objects as different filters have different OIDs
    # oids = tablePD['oid'].unique()
    # if len(oids) > 1:
    #     print('Multiple objects found. First object used.')
    #     tablePD = tablePD[tablePD['oid'] == oids[0]]

    fi = pd.Series({"zg": "g", "zr": "r", "zi": "i"})  # map filter ID from ZTF code (used as key in output)
    tablePD['filterID'] = tablePD['filtercode'].map(fi)

    return status, tablePD


def refZTFobj(coordinates: CoordClass, mjd, filefracday, mag):
    magLow = mag - config.magWindow
    magHigh = mag + config.magWindow
    timeLow = mjd - config.timeWindow
    timeHigh = mjd + config.timeWindow
    response = getLightCurveDataZTF(coordinates, config.regionSize,
                                    column_filters=(f'TIME={timeLow}{delim}{timeHigh}',
                                                    f'MAG={magLow}{delim}{magHigh}'))

    if response[0]:
        df = response[1]
        df = df[df['filefracday'] == filefracday]
        return df if df is not None else False
    pass


def getOIDZTFinfo(oid: str):
    # response = []
    # oidListFile = OIDListFile()
    # for oid in oidList:
    inList = oidListFile.inList(oid)
    if len(inList):
        response = inList
    else:
        oidSQL = 'oid=' + str(oid)
        query = f"SELECT * FROM ztf_objects_dr5 WHERE {oidSQL} ORDER BY oid"
        queryReturn = ZTFObjectQuery(query)
        if len(queryReturn):
            oidListFile.add(queryReturn)
            response = pd.DataFrame(queryReturn)
        else:
            response = None
        # response.append(response)

    return response


def getZTFOidOffsetStats(refObs):
    temp = refObs["oid"].astype(str)
    oidList = 'oid=' + temp.str.cat(sep=' OR oid=')

    query = f"SELECT * FROM ztf_objects_dr5 WHERE {oidList} ORDER BY oid"

    response = ZTFObjectQuery(query)

    final = pd.merge(pd.DataFrame(response), refObs, on="oid")
    final['offsetmag'] = final['refmag'] - refObs['mag']

    return final['offsetmag'].mean(), final['offsetmag'].std()


def ZTFObjectQuery(query):
    service = pyvo.dal.TAPService('https://irsa.ipac.caltech.edu/TAP')

    print(f'Requesting OID data from Zwicky Transient Facility. Please wait ... ', end='')
    with Spinner():
        try:
            response = service.run_async(query)
            print(f'\r{" ":66}\r ', end='')
            # print(' ', end='')
        except HTTPError as err:
            if err.code == 400:
                print('Sorry. Could not complete request.')
                return False
            else:
                return False

    return response


def getPosSDSSMedian(pos: CoordClass):
    pass


def getZTFImage(mjd, filefracday, field, filtercode, ccdid, quadrant, pos, size=None):
    # get single image
    if size is None:
        size = [5, 1, 0.1]

    hdul = []
    baseURL = 'https://irsa.ipac.caltech.edu/ibe/data/ztf/products/sci/'
    dateTime = astropy.time.Time(mjd, format='mjd')
    fracdate = str(filefracday)[-6:]
    imgtypecode = 'o'
    for i in range(len(size)):
        getURL = baseURL + dateTime.strftime('%Y/%m%d/') + fracdate + '/ztf_' + str(filefracday) \
                 + '_' + f'{field:0>6d}' + '_' + filtercode \
                 + '_' + f'c{int(ccdid, base=16):0>2d}' + '_' + imgtypecode \
                 + '_' + f'q{int(quadrant, base=16):0>1d}' + f'_sciimg.fits' \
                 + '?' + f'center={pos.ra_str()},{pos.dec_str()}&size={str(size[i])}arcmin&gzip=false'
        hdul.append(fits.open(getURL))

        # if hdul.inf != 200:  # Ensure good response is received back from IRSA
        #     return config.badResponse

        # config.LClog.debug(hdul[i].info())
    return hdul, size


def getZTFImageData(mjd, filefracday, field, filtercode, ccdid, quadrant, pos, size=None):
    # get single image
    if size is None:
        size = [1]
    ztfImageProducts = ['sciimg.fits', 'mskimg.fits', 'psfcat.fits', 'sexcat.fits', 'sciimgdao.psf',
                        'sciimgdaopsfcent.fits', 'scimrefdiffimg.fits.fz', 'diffimgpsf.fits',
                        'sciimlog.txt', 'diffimlog.txt', 'log.txt']
    hdul = []
    baseURL = 'https://irsa.ipac.caltech.edu/ibe/data/ztf/products/sci/'
    dateTime = astropy.time.Time(mjd, format='mjd')
    fracdate = str(filefracday)[-6:]
    imgtypecode = 'o'
    for prod in ztfImageProducts:
        getURL = baseURL + dateTime.strftime('%Y/%m%d/') + fracdate + '/ztf_' + str(filefracday) \
                 + '_' + f'{field:0>6d}' + '_' + filtercode \
                 + '_' + f'c{int(ccdid, base=16):0>2d}' + '_' + imgtypecode \
                 + '_' + f'q{int(quadrant, base=16):0>1d}' + f'_{prod}' \
                 + '?' + f'center={pos.ra_str()},{pos.dec_str()}&size={str(size[0])}arcmin&gzip=false'
        if prod[-4:] == 'fits' or prod[-7:] == 'fits.fz':
            hdul.append(fits.open(getURL))
        elif prod[-3:] == 'psf':
            pass
            hdul.append(pd.read_fwf(getURL, widths=(14, 13, 13, 13, 13, 13)))
        elif prod[-3:] == 'txt':
            hdul.append(requests.get(getURL))

    for hdu in hdul:
        try:
            hdu.info()
        except AttributeError:
            # ignore for text data
            pass
        finally:
            pass

        # if hdul.inf != 200:  # Ensure good response is received back from IRSA
        #     return config.badResponse

        # config.LClog.debug(hdul[i].info())
    return hdul, size


def getZTFRefImage(field, filtercode, ccdid, quadrant, pos, size=None):
    # Get reference image
    if size is None:
        size = [5, 1, 0.1]
    hdul = []
    baseURL = 'https://irsa.ipac.caltech.edu/ibe/data/ztf/products/ref/'
    paddedfield = f'{field:0>6d}'
    fieldprefix = paddedfield[:3]
    paddedccdid = f'{int(ccdid, base=16):0>2d}'
    qid = f'{int(quadrant, base=16):0>1d}'

    for i in range(len(size)):
        getURL = baseURL + fieldprefix + '/field' + paddedfield + '/' + filtercode \
                 + '/ccd' + paddedccdid + '/q' + qid \
                 + '/ztf_' + paddedfield + '_' + filtercode + '_c' + paddedccdid + '_q' + qid + '_refimg' + '.fits' \
                 + '?' + f'center={pos.ra_str()},{pos.dec_str()}&size={str(size[i])}arcmin&gzip=false'
        hdul.append(fits.open(getURL))

        # config.LClog.debug(hdul[i].info())
    return hdul, size
