import sys

from src.SDSSRefs.histogram import histogram
from src.SDSSRefs.resample import resample

from src.LCExtract.LCExtract import LCExtract
from src.SDSSRefs.controller import controller


def help():
    print('Main process for Lightcurve routines\n\n'
          'Valid command line parameters are:\n'
          'lc       : To execute the lightcurve extraction application.\n'
          'sdss     : Perform SDSS reference calibration using config parameters\n'
          'resample : Parse and reset stats from the LC files for multiple filters\n'
          'hist     : Histogram data parse')


if __name__ == '__main__':

    if len(sys.argv) == 1 or sys.argv[1] == 'lc':  # by default
        LCExtract()
    else:
        execute = sys.argv[1]
        if execute == 'resample':
            resample()
        elif execute == 'sdss':
            controller()
        elif execute == 'hist':
            histogram()
        elif execute == 'h':
            help()
    exit(0)

