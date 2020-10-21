from ..util import slice_


def csrot(N, CX, INCX, CY, INCY, C, S):
    """Applies a Givens rotation to a pair of vectors x and y

    Parameters
    ----------
    N : int
        Number of elements in input vector
    CX : numpy.ndarray
        A single precision complex array, dimension (1 + (`N` - 1)*abs(`INCX`))
    INCX : int
        Storage spacing between elements of `CX`
    CY : numpy.ndarray
        A single precision complex array, dimension (1 + (`N` - 1)*abs(`INCY`))
    INCY : int
        Storage spacing between elements of `CY`
    C : numpy.single
        The Givens parameter c, with value cos(theta)
    S : numpy.single
        The Givens parameter s, with value sin(theta)

    Returns
    -------
    None

    See Also
    --------
    srot : Single-precision real Givens rotation
    crot : Single-precision complex Givens rotation
    zdrot : Double-precision complex Givens rotation

    Notes
    -----
    Online PyBLAS documentation: https://nbviewer.jupyter.org/github/timleslie/pyblas/blob/main/docs/csrot.ipynb
    Reference BLAS documentation: https://github.com/Reference-LAPACK/lapack/blob/v3.9.0/BLAS/SRC/csrot.f

    Examples
    --------
    >>> x = np.array([1+2j, 2+3j, 3+4j], dtype=np.complex64)
    >>> y = np.array([6+7j, 7+8j, 8+9j], dtype=np.complex64)
    >>> N = len(x)
    >>> incx = 1
    >>> incy = 1
    >>> theta = np.pi/2
    >>> csrot(N, x, incx, y, incy, np.cos(theta), np.sin(theta))
    >>> print(x)
    [6.+7.j 7.+8.j 8.+9.j]
    >>> print(y)
    [-1.-2.j -2.-3.j -3.-4.j]
    """
    if N <= 0:
        return
    x_slice = slice_(N, INCX)
    y_slice = slice_(N, INCY)
    X_TEMP = C * CX[x_slice] + S * CY[y_slice]
    CY[y_slice] = -S * CX[x_slice] + C * CY[y_slice]
    CX[x_slice] = X_TEMP
