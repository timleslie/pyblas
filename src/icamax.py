# > \brief \b ICAMAX
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       INTEGER FUNCTION ICAmax(N,CX,INCX)
#
#       .. Scalar Arguments ..
#       INTEGER INCX,N
#       ..
#       .. Array Arguments ..
#       COMPLEX CX(*)
#       ..
#
#
# > \par Purpose:
#  =============
# >
# > \verbatim
# >
# >    ICAMAX finds the index of the first element having maximum |Re(.)| + |Im(.)|
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
# > \param[in] CX
# > \verbatim
# >          CX is COMPLEX array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
# > \endverbatim
# >
# > \param[in] INCX
# > \verbatim
# >          INCX is INTEGER
# >         storage spacing between elements of CX
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
# >     jack dongarra, linpack, 3/11/78.
# >     modified 3/93 to return if incx <= 0.
# >     modified 12/3/93, array(1) declarations changed to array(*)
# > \endverbatim
# >
#  =====================================================================
from scabs1 import scabs1


def ICAmax(N, CX, INCX):
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
    # COMPLEX CX(*)
    #     ..
    #
    #  =====================================================================

    if N < 1 or INCX <= 0:
        return 0
    if N == 1:
        return 1

    ICAMAX = 1
    if INCX == 1:
        #
        #        code for increment equal to 1
        #
        SMAX = scabs1(CX[1])
        for I in range(1, N):
            if scabs1(CX[I]) > SMAX:
                ICAMAX = I
                SMAX = scabs1(CX[I])
    else:
        #
        # code for increment not equal to 1
        #
        IX = 1
        SMAX = scabs1(CX[1])
        IX += INCX
        for I in range(1, N):
            if scabs1(CX[IX]) > SMAX:
                ICAMAX = I
                SMAX = scabs1(CX[IX])
            IX += INCX
    return ICAMAX
