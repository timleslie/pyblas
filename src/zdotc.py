# > \brief \b ZDOTC
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       COMPLEX*16 FUNCTION ZDOTC(N,ZX,INCX,ZY,INCY)
#
#       .. Scalar Arguments ..
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
# > ZDOTC forms the dot product of two complex vectors
# >      ZDOTC = X^H * Y
# >
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
# > \param[in] ZY
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
def ZDOTC(N, ZX, INCX, ZY, INCY):
    #
    #  -- Reference BLAS level1 routine (version 3.8.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     November 2017
    #
    #     .. Scalar Arguments ..
    # INTEGER INCX,INCY,N
    #     ..
    #     .. Array Arguments ..
    # COMPLEX*16 ZX(*),ZY(*)
    #     ..
    #
    #  =====================================================================

    ZTEMP = 0
    if N <= 0:
        return 0
    if INCX == 1 and INCY == 1:
        # code for both increments equal to 1
        for I in range(N):
            ZTEMP += (ZX[I]).conjugate() * ZY[I]
    else:
        # code for unequal increments or equal increments not equal to 1
        IX = 1
        IY = 1
        if INCX < 0:
            IX = (-N + 1) * INCX + 1
        if INCY < 0:
            IY = (-N + 1) * INCY + 1
        for I in range(N):
            ZTEMP += (ZX[IX]).conjugate() * ZY[IY]
            IX += INCX
            IY += INCY
    return ZTEMP
