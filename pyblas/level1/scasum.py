import numpy as np
from ..util import slice_


def scasum(N, CX, INCX):
    """Computes the sum of absolute values of components the vector x

    Parameters
    ----------
    N : int
        Number of elements in input vector
    CX : numpy.ndarray
        A single precision complex array, dimension (1 + (`N` - 1)*abs(`INCX`))
    INCX : int
        Storage spacing between elements of `CX`

    Returns
    -------
    numpy.single

    See Also
    --------
    dzasum : Double-precision sum of absolute component values of a vector

    Notes
    -----
    Online PyBLAS documentation: https://nbviewer.jupyter.org/github/timleslie/pyblas/blob/main/docs/scasum.ipynb
    Reference BLAS documentation: https://github.com/Reference-LAPACK/lapack/blob/v3.9.0/BLAS/SRC/scasum.f

    Examples
    --------
    >>> x = np.array([1+2j, 2+3j, 3+4j], dtype=np.complex64)
    >>> N = len(x)
    >>> incx = 1
    >>> print(scasum(N, x, incx)
    15.0
    """
    if N <= 0:
        return 0
    s = slice_(N, INCX)
    return (np.abs(CX[s].real) + np.abs(CX[s].imag)).sum()
