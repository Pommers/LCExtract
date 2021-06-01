import itertools
import math
import sys
import threading
import time
from numbers import Number

from LCExtract import config


class Namespace(object):
    """
    helps referencing object in a dictionary as dict.key instead of dict['key']
    """
    def __init__(self, adict):
        self.__dict__.update(adict)


class Spinner:
    busy = False
    delay = 0.2

    def __init__(self, delay=None):
        self.spinner_generator = itertools.cycle('-/|\\')
        if delay and isinstance(delay, Number):
            self.delay = delay

    def spinner_task(self):
        while self.busy:
            sys.stdout.write(next(self.spinner_generator))
            sys.stdout.flush()
            time.sleep(self.delay)
            sys.stdout.write('\b')
            sys.stdout.flush()

    def __enter__(self):
        self.busy = True
        threading.Thread(target=self.spinner_task).start()

    def __exit__(self, exception, value, tb):
        self.busy = False
        time.sleep(self.delay)
        if exception is not None:
            return False


def threeSigma(archive, filterID, mag) -> float:
    """Calculate 3 sigma confidence level from coefficients

    previously calculated coefficients used to provide a 3 sigma (y) value for a mag (x) value passed

    :param archive: archive to calculate
    :type archive:
    :param filterID: filter to calculate
    :type filterID: str
    :param mag: mag used as 'x' value
    :type mag: float
    :return: 3 sigma value
    :rtype: float
    """
    y1 = 0
    if filterID == 'g':
        c = config.medSD_c
        for o in range(len(c)):
            y1 += c[o] * mag ** (len(c) - o - 1)
    return 3 * y1
