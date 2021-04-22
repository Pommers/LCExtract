"""
Module
------
coord.py

Description
-----------
Currently a placeholder for development of any coordinate manipulation routines which are required
in extension to the Astropy coordinate module (which in fairness seems quite rich!)
"""
from astropy.coordinates import SkyCoord


def to_string(value, precision):
    return f"{value:.{precision}f}"


class CoordClass:
    def __init__(self, ra, dec):
        self.skyCoord = SkyCoord(ra, dec, frame='icrs', unit='deg')
        self.precision = 6

    def getRA(self):
        return self.skyCoord.ra.degree

    def getDEC(self):
        return self.skyCoord.dec.degree

    def ra_str(self):
        return to_string(self.skyCoord.ra.degree, self.precision)

    def dec_str(self):
        return to_string(self.skyCoord.dec.degree, self.precision)
