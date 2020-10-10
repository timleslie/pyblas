from ..util import slice_


def zcopy(N, ZX, INCX, ZY, INCY):
    """Copies a vector, x, to a vector, y.

    Parameters
    ----------
    N : int
        Number of elements in input vectors
    ZX : numpy.ndarray
        A double precision complex array, dimension (1 + (`N` - 1)*abs(`INCX`))
    INCX : int
        Storage spacing between elements of `ZX`
    ZY : numpy.ndarray
        A single precision complex array, dimension (1 + (`N` - 1)*abs(`INCY`))
    INCY : double
        Storage spacing between elements of `ZY`

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
    dcopy : Double-precision real copy
    ccopy : Single-precision complex copy

    Notes
    -----
    Reference BLAS documentation: https://github.com/Reference-LAPACK/lapack/blob/v3.9.0/BLAS/SRC/zcopy.f

    Examples
    --------
    >>> x = np.array([1+2j, 2+3j, 3+4j], dtype=np.complex128)
    >>> y = np.array([6+7j, 7+8j, 8+9j], dtype=np.complex128)
    >>> N = len(x)
    >>> incx = 1
    >>> incy = 1
    >>> zcopy(N, x, incx, y, incy)
    >>> print(y)
    [1.+2.j 2.+3.j 3.+4.j]
    """

    if N <= 0:
        return
    ZY[slice_(N, INCY)] = ZX[slice_(N, INCX)]
