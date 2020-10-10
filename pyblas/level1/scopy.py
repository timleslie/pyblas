from ..util import slice_


def scopy(N, SX, INCX, SY, INCY):
    """Copies a vector, x, to a vector, y.

    Parameters
    ----------
    N : int
        Number of elements in input vectors
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

    Raises
    ------
    Exception
        If N <= 0

    See Also
    --------
    dcopy : Double-precision real copy
    ccopy : Single-precision complex copy
    zcopy : Double-precision complex copy

    Notes
    -----
    Reference BLAS documentation: https://github.com/Reference-LAPACK/lapack/blob/v3.9.0/BLAS/SRC/scopy.f

    Examples
    --------
    >>> x = np.array([1, 2, 3], dtype=np.single)
    >>> y = np.array([6, 7, 8], dtype=np.single)
    >>> N = len(x)
    >>> incx = 1
    >>> incy = 1
    >>> scopy(N, x, incx, y, incy)
    >>> print(y)
    [1. 2. 3.]
    """

    if N <= 0:
        return
    SY[slice_(N, INCY)] = SX[slice_(N, INCX)]
