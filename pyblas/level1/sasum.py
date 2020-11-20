from ..util import slice_


def sasum(N, SX, INCX):
    """Computes the sum of absolute values of elements of the vector x

    Parameters
    ----------
    N : int
        Number of elements in input vector
    SX : numpy.ndarray
        A single precision real array, dimension (1 + (`N` - 1)*abs(`INCX`))
    INCX : int
        Storage spacing between elements of `SX`

    Returns
    -------
    numpy.single

    See Also
    --------
    dasum : Double-precision sum of absolute values

    Notes
    -----
    Online PyBLAS documentation: https://nbviewer.jupyter.org/github/timleslie/pyblas/blob/main/docs/sasum.ipynb
    Reference BLAS documentation: https://github.com/Reference-LAPACK/lapack/blob/v3.9.0/BLAS/SRC/sasum.f

    Examples
    --------
    >>> x = np.array([1, 2, 3], dtype=np.single)
    >>> N = len(x)
    >>> incx = 1
    >>> print(sasum(N, x, incx)
    6.
    """
    if N <= 0:
        return 0
    return SX[slice_(N, INCX)].sum()
