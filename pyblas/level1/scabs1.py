import numpy as np


def scabs1(Z):
    """Computes the sum of the absolute value of real and imaginary components of z

    Parameters
    ----------
    Z : numpy.complex64
        The single-precision complex number z

    Returns
    -------
    numpy.single

    See Also
    --------
    dcabs1 : Double-precision absolute sum of components

    Notes
    -----
    Online PyBLAS documentation: https://nbviewer.jupyter.org/github/timleslie/pyblas/blob/main/docs/scabs1.ipynb
    Reference BLAS documentation: https://github.com/Reference-LAPACK/lapack/blob/v3.9.0/BLAS/SRC/scabs1.f

    Examples
    --------
    >>> z = np.complex64(-1 + 2j)
    >>> print(scabs1(Z)
    3.0
    """
    return np.abs(Z.real) + np.abs(Z.imag)
