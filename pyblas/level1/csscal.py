from ..util import slice_


def csscal(N, SA, CX, INCX):
    """Scales a vector, x, by a constant alpha

    Parameters
    ----------
    N : int
        Number of elements in input vector
    SA : numpy.single
        Specifies the scalar alpha
    CX : numpy.ndarray
        A single precision complex array, dimension (1 + (`N` - 1)*abs(`INCX`))
    INCX : int
        Storage spacing between elements of `CX`

    Returns
    -------
    None

    See Also
    --------
    sscal : Single-precision real scaling by a real constant
    dscal : Double-precision real scaling by a real constant
    cscal : Single-precision complex scaling by a complex constant
    zscal : Double-precision complex scaling by a complex constant
    zdscal : Double-precision complex scaling by a real constant

    Notes
    -----
    Online PyBLAS documentation: https://nbviewer.jupyter.org/github/timleslie/pyblas/blob/main/docs/csscal.ipynb
    Reference BLAS documentation: https://github.com/Reference-LAPACK/lapack/blob/v3.9.0/BLAS/SRC/csscal.f

    Examples
    --------
    >>> x = np.array([1+2j, 2+3j, 3+4j], dtype=np.complex64)
    >>> N = len(x)
    >>> alpha = 5
    >>> incx = 1
    >>> csscal(N, alpha, x, incx)
    >>> print(x)
    [ 5.+10.j 10.+15.j 15.+20.j]
    """
    if N <= 0:
        return
    CX[slice_(N, INCX)] *= SA
