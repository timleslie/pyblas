import numpy as np
from ..util import slice_


def dznrm2(N, X, INCX):
    """Computes the Euclidean norm of the vector x

    Parameters
    ----------
    N : int
        Number of elements in input vector
    X : numpy.ndarray
        A double precision complex array, dimension (1 + (`N` - 1)*abs(`INCX`))
    INCX : int
        Storage spacing between elements of `X`

    Returns
    -------
    numpy.double

    See Also
    --------
    snrm2 : Single-precision real euclidean norm
    dnrm2 : Double-precision real euclidean norm
    scnrm2 : Real-precision complex euclidean norm

    Notes
    -----
    Online PyBLAS documentation: https://nbviewer.jupyter.org/github/timleslie/pyblas/blob/main/docs/dznrm2.ipynb
    Reference BLAS documentation: https://github.com/Reference-LAPACK/lapack/blob/v3.9.0/BLAS/SRC/dznrm2.f

    Examples
    --------
    >>> x = np.array([1+2j, 2+3j, 3+4j], dtype=np.complex64)
    >>> N = len(x)
    >>> incx = 1
    >>> print(dznrm2(N, x, incx)
    6.557438524302
    """
    if N <= 0:
        return 0
    return np.sqrt((X[slice_(N, INCX)].conj() * X[slice_(N, INCX)]).sum().real)
