from ..util import slice_


def caxpy(N, CA, CX, INCX, CY, INCY):
    """Adds a vector x times a constant alpha to a vector y

    Parameters
    ----------
    N : int
        Number of elements in input vector
    CA : numpy.complex64
        Specifies the scalar alpha
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
    saxpy : Single-precision real adding a scaled vector to a vector
    daxpy : Double-precision real adding a scaled vector to a vector
    zaxpy : Double-precision complex adding a scaled vector to a vector

    Notes
    -----
    Online PyBLAS documentation: https://nbviewer.jupyter.org/github/timleslie/pyblas/blob/main/docs/caxpy.ipynb
    Reference BLAS documentation: https://github.com/Reference-LAPACK/lapack/blob/v3.9.0/BLAS/SRC/caxpy.f

    Examples
    --------
    >>> x = np.array([1, 2, 3], dtype=np.complex64)
    >>> y = np.array([6, 7, 8], dtype=np.complex64)
    >>> N = len(x)
    >>> alpha = 5j
    >>> incx = 1
    >>> caxpy(N, alpha, x, incx, y, incy)
    >>> print(y)
    [6.+5.j 7.+10.j 8.+15.j]
    """
    if N <= 0:
        return
    CY[slice_(N, INCY)] += CA * CX[slice_(N, INCX)]
