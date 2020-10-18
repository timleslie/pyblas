from ..util import slice_


def cswap(N, CX, INCX, CY, INCY):
    """Swaps the contents of a vector x with a vector y

    Parameters
    ----------
    N : int
        Number of elements in input vector
    CX : numpy.ndarray
        A single precision complex array, dimension (1 + (`N` - 1)*abs(`INCX`))
    INCX : int
        Storage spacing between elements of `CX`
    CY : numpy.ndarray
        A single precision complex array, dimension (1 + (`N` - 1)*abs(`INCY`))
    INCY : int
        Storage spacing between elements of `CY`

    Returns
    -------
    None

    See Also
    --------
    sswap : Single-precision real swap two vectors
    dswap : Double-precision real swap two vectors
    zswap : Double-precision complex swap two vectors

    Notes
    -----
    Online PyBLAS documentation: https://nbviewer.jupyter.org/github/timleslie/pyblas/blob/main/docs/cswap.ipynb
    Reference BLAS documentation: https://github.com/Reference-LAPACK/lapack/blob/v3.9.0/BLAS/SRC/cswap.f

    Examples
    --------
    >>> x = np.array([1+2j, 2+3j, 3+4j], dtype=np.complex64)
    >>> y = np.array([6+7j, 7+8j, 8+9j], dtype=np.complex64)
    >>> N = len(x)
    >>> incx = 1
    >>> cswap(N, x, incx, y, incy)
    >>> print(x)
    [6+7j, 7+8j, 8+9j]
    >>> print(y)
    [1+2j, 2+3j, 3+4j]
    """
    if N <= 0:
        return
    x_slice = slice_(N, INCX)
    y_slice = slice_(N, INCY)
    X_TEMP = CX[x_slice]
    CX[x_slice] = CY[y_slice]
    CY[x_slice] = X_TEMP
