# > \brief \b DSDOT
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def DSDOT(N,SX,INCX,SY,INCY)
#
#       .. Scalar Arguments ..
#       INTEGER INCX,INCY,N
#       ..
#       .. Array Arguments ..
#       REAL SX(*),SY(*)
#       ..
#
#    AUTHORS
#    =======
#    Lawson, C. L., (JPL), Hanson, R. J., (SNLA),
#    Kincaid, D. R., (U. of Texas), Krogh, F. T., (JPL)
#
#
# > \par Purpose:
#  =============
# >
# > \verbatim
# >
# > Compute the inner product of two vectors with extended
# > precision accumulation and result.
# >
# > Returns D.P. dot product accumulated in D.P., for S.P. SX and SY
# > DSDOT = sum for I = 0 to N-1 of  SX(LX+I*INCX) * SY(LY+I*INCY),
# > where LX = 1 if INCX >= 0, else LX = 1+(1-N)*INCX, and LY is
# > defined in a similar way using INCY.
# > \endverbatim
#
#  Arguments:
#  ==========
#
# > \param[in] N
# > \verbatim
# >          N is INTEGER
# >         number of elements in input vector(s)
# > \endverbatim
# >
# > \param[in] SX
# > \verbatim
# >          SX is REAL array, dimension(N)
# >         single precision vector with N elements
# > \endverbatim
# >
# > \param[in] INCX
# > \verbatim
# >          INCX is INTEGER
# >          storage spacing between elements of SX
# > \endverbatim
# >
# > \param[in] SY
# > \verbatim
# >          SY is REAL array, dimension(N)
# >         single precision vector with N elements
# > \endverbatim
# >
# > \param[in] INCY
# > \verbatim
# >          INCY is INTEGER
# >         storage spacing between elements of SY
# > \endverbatim
# >
# > \result DSDOT
# > \verbatim
# >          DSDOT is DOUBLE PRECISION
# >         DSDOT  double precision dot product (zero if N<=0)
# > \endverbatim
#
#  Authors:
#  ========
#
# > \author Univ. of Tennessee
# > \author Univ. of California Berkeley
# > \author Univ. of Colorado Denver
# > \author NAG Ltd.
#
# > \date December 2016
#
# > \ingroup double_blas_level1
#
# > \par Further Details:
#  =====================
# >
# > \verbatim
# > \endverbatim
#
# > \par References:
#  ================
# >
# > \verbatim
# >
# >
# >  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
# >  Krogh, Basic linear algebra subprograms for Fortran
# >  usage, Algorithm No. 539, Transactions on Mathematical
# >  Software 5, 3 (September 1979), pp. 308-323.
# >
# >  REVISION HISTORY  (YYMMDD)
# >
# >  791001  DATE WRITTEN
# >  890831  Modified array declarations.  (WRB)
# >  890831  REVISION DATE from Version 3.2
# >  891214  Prologue converted to Version 4.0 format.  (BAB)
# >  920310  Corrected definition of LX in DESCRIPTION.  (WRB)
# >  920501  Reformatted the REFERENCES section.  (WRB)
# >  070118  Reformat to LAPACK style (JL)
# > \endverbatim
# >
#  =====================================================================
import numpy as np
from ..util import slice_


def dsdot(N, SX, INCX, SY, INCY):
    """Computes the dot-product of a vector x and a vector y.

    Parameters
    ----------
    N : int
        Number of elements in input vectors
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
    numpy.double

    See Also
    --------
    sdot : Single-precision real dot product
    sdsdot : Single-precision real dot product (computed in double precision, returned as single precision)
    ddot : Double-precision real dot product
    cdotu : Single-precision complex dot product
    cdotc : Single-precision complex conjugate dot product
    zdotu : Double-precision complex dot product
    zdotc : Double-precision complex conjugate dot product

    Notes
    -----
    Online PyBLAS documentation: https://nbviewer.jupyter.org/github/timleslie/pyblas/blob/main/docs/dsdot.ipynb
    Reference BLAS documentation: https://github.com/Reference-LAPACK/lapack/blob/v3.9.0/BLAS/SRC/dsdot.f

    Examples
    --------
    >>> x = np.array([1, 2, 3], dtype=np.single)
    >>> y = np.array([6, 7, 8], dtype=np.single)
    >>> N = len(x)
    >>> incx = 1
    >>> incy = 1
    >>> dsdot(N, x, incx, y, incy)
    44.0
    """
    if N <= 0:
        return 0
    return (
        SX[slice_(N, INCX)].astype(np.double) * SY[slice_(N, INCY)].astype(np.double)
    ).sum()
