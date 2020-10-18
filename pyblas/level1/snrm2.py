import numpy as np
from ..util import slice_


def snrm2(N, X, INCX):
    """Computes the Euclidean norm of the vector x

    Parameters
    ----------
    N : int
        Number of elements in input vector
    X : numpy.ndarray
        A single precision real array, dimension (1 + (`N` - 1)*abs(`INCX`))
    INCX : int
        Storage spacing between elements of `X`

    Returns
    -------
    numpy.single

    See Also
    --------
    dnrm2 : Double-precision real euclidean norm
    scnrm2 : Single-precision complex euclidean norm
    zdnrm2 : Double-precision complex euclidean norm

    Notes
    -----
    Online PyBLAS documentation: https://nbviewer.jupyter.org/github/timleslie/pyblas/blob/main/docs/snrm2.ipynb
    Reference BLAS documentation: https://github.com/Reference-LAPACK/lapack/blob/v3.9.0/BLAS/SRC/snrm2.f

    Examples
    --------
    >>> x = np.array([1, 2, 3], dtype=np.single)
    >>> N = len(x)
    >>> incx = 1
    >>> print(snrm2(N, x, incx)
    3.7416575
    """
    if N <= 0:
        return 0
    # Note: This implementaiton suffers from potential overflow errors for large vector values.
    # More sophisticated implementations can avoid this with appropriate scaling applied before
    # taking the square of large values.
    return np.sqrt((X[slice_(N, INCX)] ** 2).sum())
