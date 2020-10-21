from ..util import slice_


def drot(N, DX, INCX, DY, INCY, C, S):
    """Applies a Givens rotation to a pair of vectors x and y

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
    C : numpy.double
        The Givens parameter c, with value cos(theta)
    S : numpy.double
        The Givens parameter s, with value sin(theta)

    Returns
    -------
    None

    See Also
    --------
    srot : Single-precision real Givens rotation
    csrot : Single-precision complex Givens rotation
    zdrot : Double-precision complex Givens rotation

    Notes
    -----
    Online PyBLAS documentation: https://nbviewer.jupyter.org/github/timleslie/pyblas/blob/main/docs/drot.ipynb
    Reference BLAS documentation: https://github.com/Reference-LAPACK/lapack/blob/v3.9.0/BLAS/SRC/drot.f

    Examples
    --------
    >>> x = np.array([1, 2, 3], dtype=np.double)
    >>> y = np.array([6, 7, 8], dtype=np.double)
    >>> N = len(x)
    >>> incx = 1
    >>> incy = 1
    >>> theta = np.pi/2
    >>> drot(N, x, incx, y, incy, np.cos(theta), np.sin(theta))
    >>> print(x)
    [6. 7. 8.]
    >>> print(y)
    [-1. -2. -3.]
    """
    if N <= 0:
        return
    x_slice = slice_(N, INCX)
    y_slice = slice_(N, INCY)
    X_TEMP = C * DX[x_slice] + S * DY[y_slice]
    DY[y_slice] = -S * DX[x_slice] + C * DY[y_slice]
    DX[x_slice] = X_TEMP
