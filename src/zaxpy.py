# > \brief \b ZAXPY
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def ZAXPY(N,ZA,ZX,INCX,ZY,INCY)
#
#       .. Scalar Arguments ..
#       COMPLEX*16 ZA
#       INTEGER INCX,INCY,N
#       ..
#       .. Array Arguments ..
#       COMPLEX*16 ZX(*),ZY(*)
#       ..
#
#
# > \par Purpose:
#  =============
# >
# > \verbatim
# >
# >    ZAXPY constant times a vector plus a vector.
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
# > \param[in] ZA
# > \verbatim
# >          ZA is COMPLEX*16
# >           On entry, ZA specifies the scalar alpha.
# > \endverbatim
# >
# > \param[in] ZX
# > \verbatim
# >          ZX is COMPLEX*16 array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
# > \endverbatim
# >
# > \param[in] INCX
# > \verbatim
# >          INCX is INTEGER
# >         storage spacing between elements of ZX
# > \endverbatim
# >
# > \param[in,out] ZY
# > \verbatim
# >          ZY is COMPLEX*16 array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
# > \endverbatim
# >
# > \param[in] INCY
# > \verbatim
# >          INCY is INTEGER
# >         storage spacing between elements of ZY
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
# > \ingroup complex16_blas_level1
#
# > \par Further Details:
#  =====================
# >
# > \verbatim
# >
# >     jack dongarra, 3/11/78.
# >     modified 12/3/93, array(1) declarations changed to array(*)
# > \endverbatim
# >
#  =====================================================================
from dcabs1 import dcabs1


def ZAXPY(N, ZA, ZX, INCX, ZY, INCY):
    #
    #  -- Reference BLAS level1 routine (version 3.8.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     November 2017
    #
    #     .. Scalar Arguments ..
    # COMPLEX*16 ZA
    # INTEGER INCX,INCY,N
    #     ..
    #     .. Array Arguments ..
    # COMPLEX*16 ZX(*),ZY(*)
    #     ..
    #
    #  =====================================================================

    if N <= 0:
        return
    if dcabs1(ZA) == 0:
        return
    if INCX == 1 and INCY == 1:
        # code for both increments equal to 1
        for I in range(N):
            ZY[I] += ZA * ZX[I]
    else:
        # code for unequal increments or equal increments not equal to 1
        IX = 1
        IY = 1
        if INCX < 0:
            IX = (-N + 1) * INCX + 1
        if INCY < 0:
            IY = (-N + 1) * INCY + 1
        for I in range(N):
            ZY[IY] += ZA * ZX[IX]
            IX += INCX
            IY += INCY
