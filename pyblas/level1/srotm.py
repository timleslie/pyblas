from ..util import slice_


def srotm(N, SX, INCX, SY, INCY, SPARAM):
    """Applies a modified Givens rotation to a pair of vectors x and y

    Parameters
    ----------
    N : int
        Number of elements in input vector
    SX : numpy.ndarray
        A single precision real array, dimension (1 + (`N` - 1)*abs(`INCX`))
    INCX : int
        Storage spacing between elements of `SX`
    SY : numpy.ndarray
        A single precision real array, dimension (1 + (`N` - 1)*abs(`INCY`))
    INCY : int
        Storage spacing between elements of `SY`
    SPARAM : list of numpy.single
        SPARAM[0]: A flag indicating the form of the rotation matrix H
            -2: The identity matrix
            -1: A full matrix
             0: A matrix with h_11 = h_22 = 1
             1: A matrix with h_12 = h_21 = 1
        SPARAM[1-4]: The values of h_11, h_12, h_21, and h_22

    Returns
    -------
    None

    See Also
    --------
    drotm: Double-precision real modified Givens rotation

    Notes
    -----
    Online PyBLAS documentation: https://nbviewer.jupyter.org/github/timleslie/pyblas/blob/main/docs/srotm.ipynb
    Reference BLAS documentation: https://github.com/Reference-LAPACK/lapack/blob/v3.9.0/BLAS/SRC/srotm.f

    Examples
    --------
    >>> x = np.array([1, 2, 3], dtype=np.single)
    >>> y = np.array([6, 7, 8], dtype=np.single)
    >>> N = len(x)
    >>> incx = 1
    >>> incy = 1
    >>> theta = np.pi/2
    >>> srot(N, x, incx, y, incy, np.cos(theta), np.sin(theta))
    >>> print(x)
    [6. 7. 8.]
    >>> print(y)
    [-1. -2. -3.]
    """
    [SFLAG, SH11, SH21, SH12, SH22] = SPARAM
    if N <= 0:
        return
    if SFLAG == -2:  # Identity matrix
        return

    sx = slice_(N, INCX)
    sy = slice_(N, INCY)
    if SFLAG == -1:  # Full matrix
        TMP_X = SH11 * SX[sx] + SH12 * SY[sy]
        SY[sy] = SH21 * SX[sx] + SH22 * SY[sy]
        SX[sx] = TMP_X
    elif SFLAG == 0:  # Off-diagonal
        TMP_X = SX[sx] + SH12 * SY[sy]
        SY[sy] = SH21 * SX[sx] + SY[sy]
        SX[sx] = TMP_X
    elif SFLAG == 1:  # Diagonal
        TMP_X = SH11 * SX[sx] + SY[sy]
        SY[sy] = -SX[sx] + SH22 * SY[sy]
        SX[sx] = TMP_X
