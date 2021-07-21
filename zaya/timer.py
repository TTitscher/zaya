from time import perf_counter
from tabulate import tabulate

_timings = []


def list_timings(printer=None):
    printer = print if printer is None else printer
    out = tabulate(_timings, headers=["what", "time [ms]"])
    printer(out)


class TTimer:
    def __init__(self, what=""):
        self.what = what

    def __enter__(self):
        self.start = perf_counter()
        return self

    def __exit__(self, *args):
        ms = (perf_counter() - self.start) * 1000
        _timings.append((self.what, ms))
