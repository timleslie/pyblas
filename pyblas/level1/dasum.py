from ..util import slice_


def dasum(N, DX, INCX):
    """Computes the sum of absolute values of elements of the vector x

    Parameters
    ----------
    N : int
        Number of elements in input vector
    DX : numpy.ndarray
        A double precision real array, dimension (1 + (`N` - 1)*abs(`INCX`))
    INCX : int
        Storage spacing between elements of `DX`

    Returns
    -------
    numpy.single

    See Also
    --------
    dasum : Double-precision sum of absolute values

    Notes
    -----
    Online PyBLAS documentation: https://nbviewer.jupyter.org/github/timleslie/pyblas/blob/main/docs/dasum.ipynb
    Reference BLAS documentation: https://github.com/Reference-LAPACK/lapack/blob/v3.9.0/BLAS/SRC/dasum.f

    Examples
    --------
    >>> x = np.array([1, 2, 3], dtype=np.double)
    >>> N = len(x)
    >>> incx = 1
    >>> print(dasum(N, x, incx)
    6.
    """
    if N <= 0:
        return 0
    return DX[slice_(N, INCX)]
