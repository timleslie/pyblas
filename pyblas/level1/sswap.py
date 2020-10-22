from ..util import slice_


def sswap(N, SX, INCX, SY, INCY):
    """Swaps the contents of a vector x with a vector y

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

    Returns
    -------
    None

    See Also
    --------
    dswap : Double-precision real swap two vectors
    cswap : Single-precision complex swap two vectors
    zswap : Double-precision complex swap two vectors

    Notes
    -----
    Online PyBLAS documentation: https://nbviewer.jupyter.org/github/timleslie/pyblas/blob/main/docs/sswap.ipynb
    Reference BLAS documentation: https://github.com/Reference-LAPACK/lapack/blob/v3.9.0/BLAS/SRC/sswap.f

    Examples
    --------
    >>> x = np.array([1, 2, 3], dtype=np.single)
    >>> y = np.array([6, 7, 8], dtype=np.single)
    >>> N = len(x)
    >>> incx = 1
    >>> sswap(N, x, incx, y, incy)
    >>> print(x)
    [6., 7., 8.]
    >>> print(y)
    [1., 2., 3.]
    """
    if N <= 0:
        return
    x_slice = slice_(N, INCX)
    y_slice = slice_(N, INCY)
    X_TEMP = SX[x_slice].copy()
    SX[x_slice] = SY[y_slice]
    SY[x_slice] = X_TEMP
