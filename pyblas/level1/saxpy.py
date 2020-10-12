from ..util import slice_


def saxpy(N, SA, SX, INCX, SY, INCY):
    """Adds a vector x times a constant alpha to a vector y

    Parameters
    ----------
    N : int
        Number of elements in input vector
    SA : numpy.single
        Specifies the scalar alpha
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
    daxpy : Double-precision real adding a scaled vector to a vector
    caxpy : Single-precision complex adding a scaled vector to a vector
    zaxpy : Double-precision complex adding a scaled vector to a vector

    Notes
    -----
    Online PyBLAS documentation: https://nbviewer.jupyter.org/github/timleslie/pyblas/blob/main/docs/saxpy.ipynb
    Reference BLAS documentation: https://github.com/Reference-LAPACK/lapack/blob/v3.9.0/BLAS/SRC/saxpy.f

    Examples
    --------
    >>> x = np.array([1, 2, 3], dtype=np.single)
    >>> y = np.array([6, 7, 8], dtype=np.single)
    >>> N = len(x)
    >>> alpha = 5
    >>> incx = 1
    >>> saxpy(N, alpha, x, incx, y, incy)
    >>> print(y)
    [11. 17. 23.]
    """
    if N <= 0:
        return
    SY[slice_(N, INCY)] += SA * SX[slice_(N, INCX)]
