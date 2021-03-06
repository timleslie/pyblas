from ..util import slice_


def zdotu(N, ZX, INCX, ZY, INCY):
    """Computes the dot-product of a vector x and a vector y.

    Parameters
    ----------
    N : int
        Number of elements in input vectors
    ZX : numpy.ndarray
        A double precision complex array, dimension (1 + (`N` - 1)*abs(`INCX`))
    INCX : int
        Storage spacing between elements of `ZX`
    ZY : numpy.ndarray
        A double precision complex array, dimension (1 + (`N` - 1)*abs(`INCY`))
    INCY : int
        Storage spacing between elements of `ZY`

    Returns
    -------
    numpy.double

    See Also
    --------
    sdot : Single-precision real dot product
    dsdot : Single-precision real dot product (computed in double precision, returned as double precision)
    sdsdot : Single-precision real dot product (computed in double precision, returned as single precision)
    ddot : Double-precision real dot product
    cdotu : Single-precision complex dot product
    cdotc : Single-precision complex conjugate dot product
    zdotc : Double-precision complex conjugate dot product

    Notes
    -----
    Online PyBLAS documentation: https://nbviewer.jupyter.org/github/timleslie/pyblas/blob/main/docs/zdotu.ipynb
    Reference BLAS documentation: https://github.com/Reference-LAPACK/lapack/blob/v3.9.0/BLAS/SRC/zdotu.f

    Examples
    --------
    >>> x = np.array([1+2j, 2+3j, 3+4j], dtype=np.complex128)
    >>> y = np.array([6+7j, 7+8j, 8+9j], dtype=np.complex128)
    >>> N = len(x)
    >>> incx = 1
    >>> incy = 1
    >>> zdotu(N, x, incx, y, incy)
    (-30+115j)
    """
    if N <= 0:
        return 0
    return (ZX[slice_(N, INCX)] * ZY[slice_(N, INCY)]).sum()
