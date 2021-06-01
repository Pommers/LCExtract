"""outputStats.py

Module
------
outputStats.py - for generating figure of SDSS / ZTF reference plots

Description
-----------
Setup 3 plots 1) ZTF scatter 2)SDSS Total objects 3) ZTF mean detections per object


"""
import os

import numpy as np
from astropy.io import ascii
from astropy.table import Table
from matplotlib import pyplot as plt

from SDSSRefs import config


def bestFit(x, y, ax):
    a = config.minMag
    b = config.maxMag

    steps = 100
    x1 = np.linspace(a, b, steps)

    deg = 8
    c = np.polyfit(x, y, deg, rcond=None, full=False, w=None, cov=False)
    print(c)
    y1 = 0
    for o in range(deg+1):
        y1 += c[o] * x1 ** (deg - o)

    ax.plot(x1, y1, "-")


class Stats:
    def __init__(self):
        statsCols = ('mag',
                     'totalSamplesMag',
                     'gMedianOfSD',
                     'gSDofSD',
                     'gMeanOfSD',
                     'gUsedSamplesMag',
                     'gMeanSamplesLC',
                     'rMedianOfSD',
                     'rSDofSD',
                     'rMeanOfSD',
                     'rUsedSamplesMag',
                     'rMeanSamplesLC'
                     )
        self.stats = Table(names=statsCols)

    def plot(self):
        fig, ax = plt.subplots(3, 1, sharex=True, figsize=(8, 11), dpi=300)

        for f in 'gr':
            # plot Scatter of ZTF lightcurve samples based on SD of SD
            ax[0].plot(self.stats['mag'], self.stats[f'{f}MedianOfSD'],
                       label=f'Median of {f}-band', c=f, marker='.')
            # ax[0].errorbar(self.stats['mag'], self.stats['medianOfSD'], yerr=self.stats['SDofSD'],
            #                c='gray', linestyle='None')
            upBound = self.stats[f'{f}MedianOfSD'] + self.stats[f'{f}SDofSD']
            lowBound = self.stats[f'{f}MedianOfSD'] - self.stats[f'{f}SDofSD']
            ax[0].fill_between(self.stats['mag'], upBound, lowBound, color=f, alpha=0.2,
                               label=f'{f}-band sample {chr(963)}')
            ax[0].set_ylabel(f'Std Dev ({chr(963)})', fontsize=10)
            # ax[0].semilogy()
            ax[0].set_title('ZTF Lightcurves Std Dev', fontsize=12)
            ax[0].legend(loc='upper center')
            ax[0].set_ylim(bottom=0)
            bestFit(self.stats['mag'], self.stats[f'{f}MedianOfSD'], ax[0])

            # plot total and used samples from SDSS
            ax[1].plot(self.stats['mag'], self.stats[f'totalSamplesMag'],
                       label=f'{f}-band Samples returned', c=f, marker='*')
            ax[1].plot(self.stats['mag'], self.stats[f'{f}UsedSamplesMag'],
                       label=f'{f}-band Samples used ({config.sigma}{chr(963)})', c=f, marker='^', alpha=0.4)
            ax[1].set_ylim(0, config.startingSampleSize * 1.1)
            ax[1].set_ylabel('Object Count', fontsize=10)
            ax[1].set_title('SDSS objects', fontsize=12)
            ax[1].legend()

            ax[2].plot(self.stats['mag'], self.stats[f'{f}MeanSamplesLC'],
                       label=f'Mean {f}-band count/object', c=f, marker='.')
            ax[2].set_xlabel('Mag', fontsize=10)
            ax[2].set_ylabel('Detection Count', fontsize=10)
            ax[2].set_title('ZTF mean detections', fontsize=12)
            ax[2].legend()

        fig.tight_layout()
        fig.show()

    def save(self):
        if not os.path.exists('data'):
            os.makedirs('data/SDSS')
        ascii.write(self.stats,
                    output=config.filename,
                    format='csv',
                    overwrite=True)

    def add_row(self, statsRow):
        self.stats.add_row(vals=statsRow)

    def removeBlanks(self):
        for f in 'gr':
            self.stats = self.stats[self.stats[f'{f}SDofSD'] != 0]

    def load(self):
        if not os.path.exists('data'):
            os.makedirs('data/SDSS')
            print('No SDSS files. Please set config \"SDSSDataRefresh\" parameter to \"True\"')
        self.stats = ascii.read(config.filename, format='csv')
