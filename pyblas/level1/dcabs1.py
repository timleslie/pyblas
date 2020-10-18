import numpy as np


def dcabs1(Z):
    """Computes the sum of the absolute value of real and imaginary components of z

    Parameters
    ----------
    Z : numpy.complex128
        The double-precision complex number z

    Returns
    -------
    numpy.double

    See Also
    --------
    scabs1 : Single-precision absolute sum of components

    Notes
    -----
    Online PyBLAS documentation: https://nbviewer.jupyter.org/github/timleslie/pyblas/blob/main/docs/dcabs1.ipynb
    Reference BLAS documentation: https://github.com/Reference-LAPACK/lapack/blob/v3.9.0/BLAS/SRC/dcabs1.f

    Examples
    --------
    >>> z = np.complex128(1 + 2j)
    >>> print(dcabs1(Z)
    3.0
    """
    return np.abs(Z.real) + np.abs(Z.imag)
