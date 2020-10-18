import numpy as np
from ..util import slice_


def sdsdot(N, SB, SX, INCX, SY, INCY):
    """Computes the dot-product of a vector x and a vector y plus a constant beta.

    Parameters
    ----------
    N : int
        Number of elements in input vectors
    SB : numpy.single
        Specifies the scalar beta
    SX : numpy.ndarray
        A single precision real array, dimension (1 + (`N` - 1)*abs(`INCX`))
    INCX : int
        Storage spacing between elements of `SX`
    SY : numpy.ndarray
        A single precision real array, dimension (1 + (`N` - 1)*abs(`INCY`))
    INCY : int
        Storage spacing between elements of `SY`

    Returns
    -------
    numpy.single

    See Also
    --------
    sdot : Single-precision real dot product
    dsdot : Single-precision real dot product (computed in double precision, returned as double precision)
    ddot : Double-precision real dot product
    cdotu : Single-precision complex dot product
    cdotc : Single-precision complex conjugate dot product
    zdotu : Double-precision complex dot product
    zdotc : Double-precision complex conjugate dot product

    Notes
    -----
    Online PyBLAS documentation: https://nbviewer.jupyter.org/github/timleslie/pyblas/blob/main/docs/sdsdot.ipynb
    Reference BLAS documentation: https://github.com/Reference-LAPACK/lapack/blob/v3.9.0/BLAS/SRC/sdsdot.f

    Examples
    --------
    >>> beta = 5
    >>> x = np.array([1, 2, 3], dtype=np.single)
    >>> y = np.array([6, 7, 8], dtype=np.single)
    >>> N = len(x)
    >>> incx = 1
    >>> incy = 1
    >>> sdot(N, beta, x, incx, y, incy)
    49.0
    """
    if N <= 0:
        return SB
    return SB + (
        SX[slice_(N, INCX)].astype(np.double) * SY[slice_(N, INCY)].astype(np.double)
    ).sum().astype(np.single)
