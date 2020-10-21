from ..util import slice_


def drotm(N, DX, INCX, DY, INCY, DPARAM):
    """Applies a modified Givens rotation to a pair of vectors x and y

    Parameters
    ----------
    N : int
        Number of elements in input vector
    DX : numpy.ndarray
        A double precision real array, dimension (1 + (`N` - 1)*abs(`INCX`))
    INCX : int
        Storage spacing between elements of `DX`
    DY : numpy.ndarray
        A double precision real array, dimension (1 + (`N` - 1)*abs(`INCY`))
    INCY : int
        Storage spacing between elements of `DY`
    DPARAM : list of numpy.double
        DPARAM[0]: A flag indicating the form of the rotation matrix H
            -2: The identity matrix
            -1: A full matrix
             0: A matrix with h_11 = h_22 = 1
             1: A matrix with h_12 = h_21 = 1
        DPARAM[1-4]: The values of h_11, h_12, h_21, and h_22

    Returns
    -------
    None

    See Also
    --------
    srotm: Single-precision real modified Givens rotation

    Notes
    -----
    Online PyBLAS documentation: https://nbviewer.jupyter.org/github/timleslie/pyblas/blob/main/docs/drotm.ipynb
    Reference BLAS documentation: https://github.com/Reference-LAPACK/lapack/blob/v3.9.0/BLAS/SRC/drotm.f

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
    [DFLAG, DH11, DH21, DH12, DH22] = DPARAM
    if N <= 0:
        return
    if DFLAG == -2:  # Identity matrix
        return

    sx = slice_(N, INCX)
    sy = slice_(N, INCY)
    if DFLAG == -1:  # Full matrix
        TMP_X = DH11 * DX[sx] + DH12 * DY[sy]
        DY[sy] = DH21 * DX[sx] + DH22 * DY[sy]
        DX[sx] = TMP_X
    elif DFLAG == 0:  # Off-diagonal
        TMP_X = DX[sx] + DH12 * DY[sy]
        DY[sy] = DH21 * DX[sx] + DY[sy]
        DX[sx] = TMP_X
    elif DFLAG == 1:  # Diagonal
        TMP_X = DH11 * DX[sx] + DY[sy]
        DY[sy] = -DX[sx] + DH22 * DY[sy]
        DX[sx] = TMP_X
