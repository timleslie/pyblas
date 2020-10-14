from ..util import slice_


def daxpy(N, DA, DX, INCX, DY, INCY):
    """Adds a vector x times a constant alpha to a vector y

    Parameters
    ----------
    N : int
        Number of elements in input vector
    DA : numpy.double
        Specifies the scalar alpha
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

    See Also
    --------
    saxpy : Single-precision real adding a scaled vector to a vector
    caxpy : Single-precision complex adding a scaled vector to a vector
    zaxpy : Double-precision complex adding a scaled vector to a vector

    Notes
    -----
    Online PyBLAS documentation: https://nbviewer.jupyter.org/github/timleslie/pyblas/blob/main/docs/daxpy.ipynb
    Reference BLAS documentation: https://github.com/Reference-LAPACK/lapack/blob/v3.9.0/BLAS/SRC/daxpy.f

    Examples
    --------
    >>> x = np.array([1, 2, 3], dtype=np.double)
    >>> y = np.array([6, 7, 8], dtype=np.double)
    >>> N = len(x)
    >>> alpha = 5
    >>> incx = 1
    >>> daxpy(N, alpha, x, incx, y, incy)
    >>> print(y)
    [11. 17. 23.]
    """
    if N <= 0:
        return
    DY[slice_(N, INCY)] += DA * DX[slice_(N, INCX)]
