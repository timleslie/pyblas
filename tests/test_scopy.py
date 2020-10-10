import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np
import numpy.testing as npt

from src.level1.scopy import scopy

dtype = np.single  # float32

def test_scopy1():
    sx = np.array([1, 2, 3], dtype=dtype)
    N = 3
    incx = 1
    sy = np.empty(3, dtype=dtype)
    incy = 1
    scopy(N, sx, incx, sy, incy)
    npt.assert_equal(sy, [1, 2, 3])

def test_scopy2():
    sx = np.array([3, -1, 2, -1, 1, -1], dtype=dtype)
    N = 3
    incx = -2
    sy = np.empty(3, dtype=dtype)
    incy = 1
    scopy(N, sx, incx, sy, incy)
    npt.assert_equal(sy, [1, 2, 3])

def test_scopy3():
    sx = np.array([1, 2, 3], dtype=dtype)
    N = 3
    incx = 1
    sy = np.zeros(5, dtype=dtype)
    incy = -2
    scopy(N, sx, incx, sy, incy)
    npt.assert_equal(sy, [3, 0, 2, 0, 1])

def test_scopy4():
    sx = np.array([3, -1, 2, -1, 1, -1], dtype=dtype)
    N = 3
    incx = -2
    sy = np.zeros(5, dtype=dtype)
    incy = -2
    scopy(N, sx, incx, sy, incy)
    npt.assert_equal(sy, [3, 0, 2, 0, 1])

def test_scopy5():
    scopy(0, None, None, None, None)
    scopy(-1, None, None, None, None)
