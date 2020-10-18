# > \brief \b SDSDOT
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def SDSDOT(N,SB,SX,INCX,SY,INCY)
#
#       .. Scalar Arguments ..
#       REAL SB
#       INTEGER INCX,INCY,N
#       ..
#       .. Array Arguments ..
#       REAL SX(*),SY(*)
#       ..
#
# > \par Purpose:
#  =============
# >
# > \verbatim
# >
# >   Compute the inner product of two vectors with extended
# >   precision accumulation.
# >
# >   Returns S.P. result with dot product accumulated in D.P.
# >   SDSDOT = SB + sum for I = 0 to N-1 of SX(LX+I*INCX)*SY(LY+I*INCY),
# >   where LX = 1 if INCX >= 0, else LX = 1+(1-N)*INCX, and LY is
# >   defined in a similar way using INCY.
# > \endverbatim
#
#  Arguments:
#  ==========
#
# > \param[in] N
# > \verbatim
# >          N is INTEGER
# >          number of elements in input vector(s)
# > \endverbatim
# >
# > \param[in] SB
# > \verbatim
# >          SB is REAL
# >          single precision scalar to be added to inner product
# > \endverbatim
# >
# > \param[in] SX
# > \verbatim
# >          SX is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
# >          single precision vector with N elements
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
# >          SY is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
# >          single precision vector with N elements
# > \endverbatim
# >
# > \param[in] INCY
# > \verbatim
# >          INCY is INTEGER
# >          storage spacing between elements of SY
# > \endverbatim
#
#  Authors:
#  ========
#
# > \author Lawson, C. L., (JPL), Hanson, R. J., (SNLA),
# > \author Kincaid, D. R., (U. of Texas), Krogh, F. T., (JPL)
#
# > \author Univ. of Tennessee
# > \author Univ. of California Berkeley
# > \author Univ. of Colorado Denver
# > \author NAG Ltd.
#
# > \date November 2017
#
# > \ingroup single_blas_level1
#
# > \par Further Details:
#  =====================
# >
# > \verbatim
# >
# >    REFERENCES
# >
# >    C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
# >    Krogh, Basic linear algebra subprograms for Fortran
# >    usage, Algorithm No. 539, Transactions on Mathematical
# >    Software 5, 3 (September 1979), pp. 308-323.
# >
# >    REVISION HISTORY  (YYMMDD)
# >
# >    791001  DATE WRITTEN
# >    890531  Changed all specific intrinsics to generic.  (WRB)
# >    890831  Modified array declarations.  (WRB)
# >    890831  REVISION DATE from Version 3.2
# >    891214  Prologue converted to Version 4.0 format.  (BAB)
# >    920310  Corrected definition of LX in DESCRIPTION.  (WRB)
# >    920501  Reformatted the REFERENCES section.  (WRB)
# >    070118  Reformat to LAPACK coding style
# > \endverbatim
# >
#  =====================================================================
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
