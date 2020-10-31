from ..util import slice_


def zdscal(N, DA, ZX, INCX):
    """Scales a vector, x, by a constant alpha

    Parameters
    ----------
    N : int
        Number of elements in input vector
    DA : numpy.double
        Specifies the scalar alpha
    ZX : numpy.ndarray
        A double precision complex array, dimension (1 + (`N` - 1)*abs(`INCX`))
    INCX : int
        Storage spacing between elements of `ZX`

    Returns
    -------
    None

    See Also
    --------
    sscal : Single-precision real scaling by a real constant
    dscal : Double-precision real scaling by a real constant
    cscal : Single-precision complex scaling by a complex constant
    csscal : Single-precision complex scaling by a real constant
    zscal : Double-precision complex scaling by a complex constant

    Notes
    -----
    Online PyBLAS documentation: https://nbviewer.jupyter.org/github/timleslie/pyblas/blob/main/docs/zdscal.ipynb
    Reference BLAS documentation: https://github.com/Reference-LAPACK/lapack/blob/v3.9.0/BLAS/SRC/zdscal.f

    Examples
    --------
    >>> x = np.array([1+2j, 2+3j, 3+4j], dtype=np.complex128)
    >>> N = len(x)
    >>> alpha = 5
    >>> incx = 1
    >>> zdscal(N, alpha, x, incx)
    >>> print(x)
    [ 5.+10.j 10.+15.j 15.+20.j]
    """
    if N <= 0:
        return
    ZX[slice_(N, INCX)] *= DA
