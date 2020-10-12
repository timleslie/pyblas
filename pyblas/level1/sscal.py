from ..util import slice_


def sscal(N, SA, SX, INCX):
    """Scales a vector, x, by a constant alpha

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

    Returns
    -------
    None

    See Also
    --------
    dscal : Double-precision real scaling by a real constant
    cscal : Single-precision complex scaling by a complex constant
    csscal : Single-precision complex scaling by a real constant
    zscal : Double-precision complex scaling by a complex constant
    zdscal : Double-precision complex scaling by a real constant

    Notes
    -----
    Online PyBLAS documentation: https://nbviewer.jupyter.org/github/timleslie/pyblas/blob/main/docs/sscal.ipynb
    Reference BLAS documentation: https://github.com/Reference-LAPACK/lapack/blob/v3.9.0/BLAS/SRC/sscal.f

    Examples
    --------
    >>> x = np.array([1, 2, 3], dtype=np.single)
    >>> N = len(x)
    >>> alpha = 5
    >>> incx = 1
    >>> sscal(N, alpha, x, incx)
    >>> print(y)
    [5. 10. 15.]
    """
    if N <= 0:
        return
    SX[slice_(N, INCX)] *= SA
