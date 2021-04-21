from unittest import TestCase

from LCExtract import config
from LCExtract.dataretrieve import getLightCurveDataPTF, AstroObjectClass


class Test(TestCase):

    def test_get_data(self):
        self.fail()

    def test_get_light_curve_data_ptf(self):
        AO = AstroObjectClass(objectName='M60-UCD1', ra=190.8998699, dec=11.5346389,
                              subtitle="Virgo cluster, Rh=24.2pc, d=16.5Mpc")
        radiusDeg = config.coneRadius
        return_type = 'VOTABLE'
        status, response = getLightCurveDataPTF(AO.pos, radiusDeg, return_type)

        assert status
