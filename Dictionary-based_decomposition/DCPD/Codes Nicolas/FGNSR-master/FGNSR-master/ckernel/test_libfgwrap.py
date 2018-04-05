import unittest
import libfgwrap as fg
import numpy as np
import numpy.testing as npt
from ctypes import c_void_p

class InitFree(unittest.TestCase):
    """Simple tests for creating / freeing data structures"""

    def testOpenClose(self):
        lenvec = 1
        weights = np.ndarray([1.0])
        verbose = False

        data_p = fg.FGNSRproj_spinit(lenvec, weights, verbose)
        self.assertIsInstance(data_p, c_void_p)
        self.assertIsNotNone(data_p.value, "Initialization failed")

        fg.FGNSRproj_spfree(data_p)
        self.assertIsNone(data_p.value, "Free'd pointer not null'd")
