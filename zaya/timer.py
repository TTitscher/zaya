from time import perf_counter


class TTimer(object):
    def __init__(self, what="", printer=None):
        self.what = what
        self.printer = print if printer is None else printer

    def __enter__(self):
        self.start = perf_counter()
        return self

    def __exit__(self, *args):
        self.s = perf_counter() - self.start
        self.ms = self.s * 1000  # millisecs
        self.printer(f"{self.what}: {self.ms:6.3f} ms")


if __name__ == "__main__":
    with TTimer("long loop"):
        for i in range(1000000):
            pass
