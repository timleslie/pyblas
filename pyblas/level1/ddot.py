from ..util import slice_


def ddot(N, DX, INCX, DY, INCY):
    """Computes the dot-product of a vector x and a vector y.

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
    numpy.double

    See Also
    --------
    sdot : Single-precision real dot product
    dsdot : Single-precision real dot product (computed in double precision, returned as double precision)
    sdsdot : Single-precision real dot product (computed in double precision, returned as single precision)
    cdotu : Single-precision complex dot product
    cdotc : Single-precision complex conjugate dot product
    zdotu : Double-precision complex dot product
    zdotc : Double-precision complex conjugate dot product

    Notes
    -----
    Online PyBLAS documentation: https://nbviewer.jupyter.org/github/timleslie/pyblas/blob/main/docs/ddot.ipynb
    Reference BLAS documentation: https://github.com/Reference-LAPACK/lapack/blob/v3.9.0/BLAS/SRC/ddot.f

    Examples
    --------
    >>> x = np.array([1, 2, 3], dtype=np.double)
    >>> y = np.array([6, 7, 8], dtype=np.double)
    >>> N = len(x)
    >>> incx = 1
    >>> incy = 1
    >>> ddot(N, x, incx, y, incy)
    44.0
    """
    if N <= 0:
        return 0
    return (DX[slice_(N, INCX)] * DY[slice_(N, INCY)]).sum()
