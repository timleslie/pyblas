from ..util import slice_


def dcopy(N, DX, INCX, DY, INCY):
    """Copies a vector, x, to a vector, y.

    Parameters
    ----------
    N : int
        Number of elements in input vectors
    DX : numpy.ndarray
        A double precision real array, dimension (1 + (`N` - 1)*abs(`INCX`))
    INCX : int
        Storage spacing between elements of `DX`
    DY : numpy.ndarray
        A double precision real array, dimension (1 + (`N` - 1)*abs(`INCY`))
    INCY : int
        Storage spacing between elements of `DY`

    Returns
    -------
    None

    Raises
    ------
    Exception
        If N <= 0

    See Also
    --------
    scopy : Single-precision real copy
    ccopy : Single-precision complex copy
    zcopy : Double-precision complex copy

    Notes
    -----
    Reference BLAS documentation: https://github.com/Reference-LAPACK/lapack/blob/v3.9.0/BLAS/SRC/dcopy.f

    Examples
    --------
    >>> x = np.array([1, 2, 3], dtype=np.double)
    >>> y = np.array([6, 7, 8], dtype=np.double)
    >>> N = len(x)
    >>> incx = 1
    >>> incy = 1
    >>> dcopy(N, x, incx, y, incy)
    >>> print(y)
    [1. 2. 3.]
    """
    if N <= 0:
        return
    DY[slice_(N, INCY)] = DX[slice_(N, INCX)]
