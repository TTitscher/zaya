import unittest
import time
import zaya

class TestPrinter:
    def __call__(self, what):
        self.what = what


class TestTimer(unittest.TestCase):
    def test_message(self):
        with zaya.TTimer("my task"):
            time.sleep(0.1)

        with zaya.TTimer("my other task"):
            pass
        tprint = TestPrinter()
        zaya.list_timings(tprint)
        zaya.list_timings()

        self.assertIn("my task", tprint.what)


if __name__ == "__main__":
    unittest.main()






