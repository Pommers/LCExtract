from unittest import TestCase

from LCExtract.ztf import refZTFobj
from LCExtract.coord import CoordClass


class Test(TestCase):
    def test_ref_ztfobj(self):
        coord = CoordClass(188.869583, -3.789194)
        response = refZTFobj(coord, 58881.5, 20200202479861, 17.57)
        assert response is not False
