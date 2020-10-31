from ..util import slice_


def cscal(N, CA, CX, INCX):
    """Scales a vector, x, by a constant alpha

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

    Returns
    -------
    None

    See Also
    --------
    sscal : Single-precision real scaling by a real constant
    dscal : Double-precision real scaling by a real constant
    csscal : Single-precision complex scaling by a real constant
    zscal : Double-precision complex scaling by a complex constant
    zdscal : Double-precision complex scaling by a real constant

    Notes
    -----
    Online PyBLAS documentation: https://nbviewer.jupyter.org/github/timleslie/pyblas/blob/main/docs/cscal.ipynb
    Reference BLAS documentation: https://github.com/Reference-LAPACK/lapack/blob/v3.9.0/BLAS/SRC/cscal.f

    Examples
    --------
    >>> x = np.array([1+2j, 2+3j, 3+4j], dtype=np.complex64)
    >>> N = len(x)
    >>> alpha = 5j
    >>> incx = 1
    >>> cscal(N, alpha, x, incx)
    >>> print(x)
    [-10. +5.j -15.+10.j -20.+15.j]
    """
    if N <= 0:
        return
    CX[slice_(N, INCX)] *= CA
