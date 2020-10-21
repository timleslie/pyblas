from ..util import slice_


def srot(N, SX, INCX, SY, INCY, C, S):
    """Applies a Givens rotation to a pair of vectors x and y

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
    C : numpy.single
        The Givens parameter c, with value cos(theta)
    S : numpy.single
        The Givens parameter s, with value sin(theta)

    Returns
    -------
    None

    See Also
    --------
    drot : Double-precision real Givens rotation
    csrot : Single-precision complex Givens rotation
    zdrot : Double-precision complex Givens rotation

    Notes
    -----
    Online PyBLAS documentation: https://nbviewer.jupyter.org/github/timleslie/pyblas/blob/main/docs/srot.ipynb
    Reference BLAS documentation: https://github.com/Reference-LAPACK/lapack/blob/v3.9.0/BLAS/SRC/srot.f

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
    if N <= 0:
        return
    x_slice = slice_(N, INCX)
    y_slice = slice_(N, INCY)
    X_TEMP = C * SX[x_slice] + S * SY[y_slice]
    SY[y_slice] = -S * SX[x_slice] + C * SY[y_slice]
    SX[x_slice] = X_TEMP
