"""SDSS.py

Module
------
SDSS.py

Description
-----------
Module for access to SDSS data. Contains SDSSdata Class

"""

from astropy.table import Table
from astroquery.sdss import SDSS

from LCExtract import config
from LCExtract.coord import CoordClass


class SDSSdata:
    def __init__(self, mag: float, pos: CoordClass):
        self.mag = mag
        self.coords = pos
        self.samples = Table
        self.filename = f'data/SDSS/SDSSg{self.mag}.csv'

    def getSamples(self):
        sampleCount = config.sampleSize
        magLow = self.mag - config.magWindow / 2
        magHigh = self.mag + config.magWindow / 2
        raHigh = self.coords.getRA() + config.regionSize
        raLow = self.coords.getRA() - config.regionSize
        decHigh = self.coords.getDEC() + config.regionSize
        decLow = self.coords.getDEC() - config.regionSize

        query = f"SELECT TOP {sampleCount} objid as Name, ra as RA, dec as DEC, g as g, r as r FROM PhotoPrimary " \
                f"WHERE ra BETWEEN {raLow} AND {raHigh} " \
                f"AND dec BETWEEN {decLow} AND {decHigh} " \
                f"AND g < 21.0 " \
                f"AND type = 6 " \
                f"AND mode = 1"
        self.samples = SDSS.query_sql(query)
        return (self.samples is not None) and (len(self.samples) > int(sampleCount/2))

