import numpy as np
from ..util import slice_


def scnrm2(N, X, INCX):
    """Computes the Euclidean norm of the vector x

    Parameters
    ----------
    N : int
        Number of elements in input vector
    X : numpy.ndarray
        A single precision complex array, dimension (1 + (`N` - 1)*abs(`INCX`))
    INCX : int
        Storage spacing between elements of `X`

    Returns
    -------
    numpy.single

    See Also
    --------
    snrm2 : Single-precision real euclidean norm
    dnrm2 : Double-precision real euclidean norm
    dznrm2 : Double-precision complex euclidean norm

    Notes
    -----
    Online PyBLAS documentation: https://nbviewer.jupyter.org/github/timleslie/pyblas/blob/main/docs/scnrm2.ipynb
    Reference BLAS documentation: https://github.com/Reference-LAPACK/lapack/blob/v3.9.0/BLAS/SRC/scnrm2.f

    Examples
    --------
    >>> x = np.array([1+2j, 2+3j, 3+4j], dtype=np.complex64)
    >>> N = len(x)
    >>> incx = 1
    >>> print(scnrm2(N, x, incx)
    6.5574384
    """
    if N <= 0:
        return 0
    # Note: This implementaiton suffers from potential overflow errors for large vector values.
    # More sophisticated implementations can avoid this with appropriate scaling applied before
    # taking the square of large values.
    return np.sqrt((X[slice_(N, INCX)].conj() * X[slice_(N, INCX)]).sum().real)
