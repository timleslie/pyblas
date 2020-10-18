from ..util import slice_


def dswap(N, DX, INCX, DY, INCY):
    """Swaps the contents of a vector x with a vector y

    Parameters
    ----------
    N : int
        Number of elements in input vector
    DX : numpy.ndarray
        A single precision real array, dimension (1 + (`N` - 1)*abs(`INCX`))
    INCX : int
        Storage spacing between elements of `DX`
    DY : numpy.ndarray
        A single precision real array, dimension (1 + (`N` - 1)*abs(`INCY`))
    INCY : int
        Storage spacing between elements of `DY`

    Returns
    -------
    None

    See Also
    --------
    sswap : Single-precision real swap two vectors
    cswap : Single-precision complex swap two vectors
    zswap : Double-precision complex swap two vectors

    Notes
    -----
    Online PyBLAS documentation: https://nbviewer.jupyter.org/github/timleslie/pyblas/blob/main/docs/dswap.ipynb
    Reference BLAS documentation: https://github.com/Reference-LAPACK/lapack/blob/v3.9.0/BLAS/SRC/dswap.f

    Examples
    --------
    >>> x = np.array([1, 2, 3], dtype=np.double)
    >>> y = np.array([6, 7, 8], dtype=np.double)
    >>> N = len(x)
    >>> incx = 1
    >>> dswap(N, x, incx, y, incy)
    >>> print(x)
    [6., 7., 8.]
    >>> print(y)
    [1., 2., 3.]
    """
    if N <= 0:
        return
    x_slice = slice_(N, INCX)
    y_slice = slice_(N, INCY)
    X_TEMP = DX[x_slice]
    DX[x_slice] = DY[y_slice]
    DY[x_slice] = X_TEMP
