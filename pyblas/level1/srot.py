# > \brief \b SROT
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def SROT(N,SX,INCX,SY,INCY,C,S)
#
#       .. Scalar Arguments ..
#       REAL C,S
#       INTEGER INCX,INCY,N
#       ..
#       .. Array Arguments ..
#       REAL SX(*),SY(*)
#       ..
#
#
# > \par Purpose:
#  =============
# >
# > \verbatim
# >
# >    applies a plane rotation.
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
# > \param[in,out] SX
# > \verbatim
# >          SX is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
# > \endverbatim
# >
# > \param[in] INCX
# > \verbatim
# >          INCX is INTEGER
# >         storage spacing between elements of SX
# > \endverbatim
# >
# > \param[in,out] SY
# > \verbatim
# >          SY is REAL array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
# > \endverbatim
# >
# > \param[in] INCY
# > \verbatim
# >          INCY is INTEGER
# >         storage spacing between elements of SY
# > \endverbatim
# >
# > \param[in] C
# > \verbatim
# >          C is REAL
# > \endverbatim
# >
# > \param[in] S
# > \verbatim
# >          S is REAL
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
# > \date November 2017
#
# > \ingroup single_blas_level1
#
# > \par Further Details:
#  =====================
# >
# > \verbatim
# >
# >     jack dongarra, linpack, 3/11/78.
# >     modified 12/3/93, array(1) declarations changed to array(*)
# > \endverbatim
# >
#  =====================================================================
from ..util import slice_


def srot(N, SX, INCX, SY, INCY, C, S):
    #
    #  -- Reference BLAS level1 routine (version 3.8.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     November 2017
    #
    #     .. Scalar Arguments ..
    #   REAL C,S
    #   INTEGER INCX,INCY,N
    #     ..
    #     .. Array Arguments ..
    #   REAL SX(*),SY(*)
    #     ..
    #
    #  =====================================================================
    #
    #     .. Local Scalars ..
    #   REAL STEMP
    #   INTEGER I,IX,IY
    #     ..
    if N <= 0:
        return
    x_slice = slice_(N, INCX)
    y_slice = slice_(N, INCY)
    X_TEMP = C * SX[x_slice] + S * SY[y_slice]
    SY[y_slice] = -S * SX[x_slice] + C * SY[y_slice]
    SX[x_slice] = X_TEMP
