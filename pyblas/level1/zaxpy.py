from ..util import slice_


def zaxpy(N, ZA, ZX, INCX, ZY, INCY):
    """Adds a vector x times a constant alpha to a vector y

    Parameters
    ----------
    N : int
        Number of elements in input vector
    ZA : numpy.complex128
        Specifies the scalar alpha
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
    None

    See Also
    --------
    saxpy : Single-precision real adding a scaled vector to a vector
    daxpy : Double-precision real adding a scaled vector to a vector
    zaxpy : Double-precision complex adding a scaled vector to a vector

    Notes
    -----
    Online PyBLAS documentation: https://nbviewer.jupyter.org/github/timleslie/pyblas/blob/main/docs/zaxpy.ipynb
    Reference BLAS documentation: https://github.com/Reference-LAPACK/lapack/blob/v3.9.0/BLAS/SRC/zaxpy.f

    Examples
    --------
    >>> x = np.array([1+2j, 2+3j, 3+4j], dtype=np.complex128)
    >>> y = np.array([6+7j, 7+8j, 8+9j], dtype=np.complex128)
    >>> N = len(x)
    >>> alpha = 5j
    >>> incx = 1
    >>> zaxpy(N, alpha, x, incx, y, incy)
    >>> print(y)
    [ -4.+12.j  -8.+18.j -12.+24.j]
    """
    if N <= 0:
        return
    ZY[slice_(N, INCY)] += ZA * ZX[slice_(N, INCX)]
