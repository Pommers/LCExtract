from pprint import pprint
from unittest import TestCase

import numpy as np
import requests
from astropy.io import ascii
from matplotlib import pyplot as plt

from LCExtract.PanSTARRS import ps1cone, ps1search, ps1metadata
from LCExtract.PanSTARRS import checklegal, mastQuery, resolve, getDetections

objName = 'M60-UCD1'


class Test(TestCase):
    def _getColumns(self):
        # set columns to return for testing purposes
        # strip blanks and weed out blank and commented-out values
        columns = """objID,raMean,decMean,nDetections,ng,nr,ni,nz,ny,
            gMeanPSFMag,rMeanPSFMag,iMeanPSFMag,zMeanPSFMag,yMeanPSFMag""".split(',')
        columns = [x.strip() for x in columns]
        columns = [x for x in columns if x and not x.startswith('#')]
        return columns

    def test_ps1cone(self):
        ra, dec = resolve(objName)
        radius = 1.0 / 3600.0  # radius = 1 arcsec
        constraints = {'nDetections.gt': 1}

        results = ps1cone(ra, dec, radius, release='dr2', columns=self._getColumns(), **constraints)
        if not results:
            self.fail('No results returned.')
        else:
            tab = ascii.read(results)
            # improve the format
            for filter in 'grizy':
                col = filter + 'MeanPSFMag'
                tab[col].format = ".4f"
                tab[col][tab[col] == -999.0] = np.nan
            pprint(tab)

    def test_ps1search(self):
        results1 = ps1search(release='dr2', columns=self._getColumns(), verbose=True, objid=str(122851876947049923))
        tab1 = ascii.read(results1)
        print(tab1)
        self.fail()

    def test_checklegal(self):
        self.fail()

    def test_ps1metadata(self):
        try:
            meta = ps1metadata('detection', 'dr2')
        except requests.exceptions.HTTPError:
            self.fail('Unable to reach server')
        except:
            self.fail()
        else:
            print(meta)
            self.assertTrue(len(meta) == 58)  # test response for number of columns returned

    def test_mast_query(self):
        self.fail()

    def test_resolve(self):
        try:
            ra, dec = resolve(objName)
        except ValueError:
            self.fail(f'Unable to resolve object {objName}')
        else:
            assert True

    def test_get_detections(self):

        filtStr = 'gri'  # need list (1, 2, 3)
        f = [i+1 for i, n in enumerate('grizy') if n in filtStr]
        print(f)

        # ra, dec = resolve(objName)
        ra, dec = 190.8998699, 11.5346389
        radius = 1.0 / 3600.0  # radius = 1 arcsec
        constraints = {'nDetections.gt': 1}

        results = ps1cone(ra, dec, radius, release='dr2', columns=self._getColumns(), **constraints)
        if not results:
            self.fail('No results returned.')
        tab = ascii.read(results)
        pprint(tab, width=400)

        # improve the format
        for filter in 'grizy':
            col = filter + 'MeanPSFMag'
            tab[col].format = ".4f"
            tab[col][tab[col] == -999.0] = np.nan
        pprint(tab, width=400)

        dTab = getDetections(tab)
        pprint(dTab, width=400)
        return dTab

    def test_lightCurvePlot(self):
        dTab = self.test_get_detections()
        # convert flux in Jy to magnitudes
        t = dTab['obsTime']
        dmag = -2.5 * np.log10(dTab['psfFlux']) + 8.90
        print(f'Total Magnitude samples (all filters): {len(dmag)}')
        dmag1 = dmag[dTab['psfQfPerfect'] > 0.90]
        print(f'Good magnitude samples (all filters): {len(dmag1)}')
        xlim = np.array([t.min(), t.max()])
        xlim = xlim + np.array([-1, 1]) * 0.02 * (xlim[1] - xlim[0])
        # flip axis direction for magnitude
        ylim = np.array([dmag1.max(), dmag1.min()])
        ylim = ylim + np.array([-1, 1]) * 0.02 * (ylim[1] - ylim[0])

        plt.rcParams.update({'font.size': 14})
        plt.figure(1, (10, 10))
        for i, filter in enumerate("grizy"):
            plt.subplot(511 + i)
            w = np.where((dTab['filter'] == filter) & (dTab['psfQfPerfect'] > 0.90))[0]
            print(f'Total Magnitude samples (filter {filter}): {len(w)}')
            plt.scatter(t[w], dmag[w], '-*')
            plt.ylabel(filter + ' [mag]')
            plt.xlim(xlim)
            plt.ylim(ylim)
            # plt.gca().invert_yaxis()
            if i == 0:
                plt.title(objName)
        plt.xlabel('Time [MJD]')
        plt.tight_layout()
        plt.show()

