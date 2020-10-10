# > \brief \b CAXPY
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def CAXPY(N,CA,CX,INCX,CY,INCY)
#
#       .. Scalar Arguments ..
#       COMPLEX CA
#       INTEGER INCX,INCY,N
#       ..
#       .. Array Arguments ..
#       COMPLEX CX(*),CY(*)
#       ..
#
#
# > \par Purpose:
#  =============
# >
# > \verbatim
# >
# >    CAXPY constant times a vector plus a vector.
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
# > \param[in] CA
# > \verbatim
# >          CA is COMPLEX
# >           On entry, CA specifies the scalar alpha.
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
# >
# > \param[in,out] CY
# > \verbatim
# >          CY is COMPLEX array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
# > \endverbatim
# >
# > \param[in] INCY
# > \verbatim
# >          INCY is INTEGER
# >         storage spacing between elements of CY
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
# > \ingroup complex_blas_level1
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
from scabs1 import scabs1


def caxpy(N, CA, CX, INCX, CY, INCY):
    #
    #  -- Reference BLAS level1 routine (version 3.8.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     November 2017
    #
    #     .. Scalar Arguments ..
    # COMPLEX CA
    # INTEGER INCX,INCY,N
    #     ..
    #     .. Array Arguments ..
    # COMPLEX CX(*),CY(*)
    #     ..
    #
    #  =====================================================================

    if N <= 0:
        return
    if INCX == 1 and INCY == 1:
        # code for both increments equal to 1
        CY[:N] += CA * CX[:N]
    else:
        # code for unequal increments or equal increments not equal to 1
        IX = 1
        IY = 1
        if INCX < 0:
            IX = (-N + 1) * INCX + 1
        if INCY < 0:
            IY = (-N + 1) * INCY + 1
        for I in range(N):
            CY[IY] = CY[IY] + CA * CX[IX]
            IX += INCX
            IY += INCY
