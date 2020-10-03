import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np
import numpy.testing as npt

from src.level1.scopy import SCOPY


def func(x):
    return x + 1


def test_scopy():
    sx = np.array([1, 2, 3])
    N = 3
    incx = 1
    sy = np.empty(3)
    incy = 1
    SCOPY(N, sx, incx, sy, incy)
    npt.assert_equal(sx, sy)
    npt.assert_equal(sx, [1, 2, 3])
    npt.assert_equal(sy, [1, 2, 3])
