import itertools
import sys
import threading
import time
from numbers import Number


class Spinner:
    busy = False
    delay = 0.1

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
