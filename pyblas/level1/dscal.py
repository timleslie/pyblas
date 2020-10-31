from ..util import slice_


def dscal(N, DA, DX, INCX):
    """Scales a vector, x, by a constant alpha

    Parameters
    ----------
    N : int
        Number of elements in input vector
    DA : numpy.single
        Specifies the scalar alpha
    DX : numpy.ndarray
        A double precision real array, dimension (1 + (`N` - 1)*abs(`INCX`))
    INCX : int
        Storage spacing between elements of `DX`

    Returns
    -------
    None

    See Also
    --------
    sscal : Single-precision real scaling by a real constant
    cscal : Single-precision complex scaling by a complex constant
    csscal : Single-precision complex scaling by a real constant
    zscal : Double-precision complex scaling by a complex constant
    zdscal : Double-precision complex scaling by a real constant

    Notes
    -----
    Online PyBLAS documentation: https://nbviewer.jupyter.org/github/timleslie/pyblas/blob/main/docs/dscal.ipynb
    Reference BLAS documentation: https://github.com/Reference-LAPACK/lapack/blob/v3.9.0/BLAS/SRC/dscal.f

    Examples
    --------
    >>> x = np.array([1, 2, 3], dtype=np.double)
    >>> N = len(x)
    >>> alpha = 5
    >>> incx = 1
    >>> dscal(N, alpha, x, incx)
    >>> print(x)
    [5. 10. 15.]
    """
    if N <= 0:
        return
    DX[slice_(N, INCX)] *= DA
