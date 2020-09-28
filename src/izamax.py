# > \brief \b IZAMAX
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       INTEGER FUNCTION IZAmax(N,ZX,INCX)
#
#       .. Scalar Arguments ..
#       INTEGER INCX,N
#       ..
#       .. Array Arguments ..
#       COMPLEX*16 ZX(*)
#       ..
#
#
# > \par Purpose:
#  =============
# >
# > \verbatim
# >
# >    IZAMAX finds the index of the first element having maximum |Re(.)| + |Im(.)|
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
# > \ingroup aux_blas
#
# > \par Further Details:
#  =====================
# >
# > \verbatim
# >
# >     jack dongarra, 1/15/85.
# >     modified 3/93 to return if incx <= 0.
# >     modified 12/3/93, array(1) declarations changed to array(*)
# > \endverbatim
# >
#  =====================================================================
from dcabs1 import dcabs1


def IZAmax(N, ZX, INCX):
    #
    #  -- Reference BLAS level1 routine (version 3.8.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     November 2017
    #
    #     .. Scalar Arguments ..
    # INTEGER INCX,N
    #     ..
    #     .. Array Arguments ..
    # COMPLEX*16 ZX(*)
    #     ..
    #
    #  =====================================================================

    if N < 1 or INCX <= 0:
        return 0
    if N == 1:
        return 1

    IZAMAX = 1
    if INCX == 1:
        #
        #        code for increment equal to 1
        #
        DMAX = dcabs1(ZX[1])
        for I in range(1, N):
            if dcabs1(ZX[I]) > DMAX:
                IZAMAX = I
                DMAX = dcabs1(ZX[I])
    else:
        #
        #        code for increment not equal to 1
        #
        IX = 1
        DMAX = dcabs1(ZX[1])
        IX += INCX
        for I in range(1, N):
            if dcabs1(ZX[IX]) > DMAX:
                IZAMAX = I
                DMAX = dcabs1(ZX[IX])
            IX += INCX
    return IZAMAX
