# > \brief \b CDOTC
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def CDOTC(N,CX,INCX,CY,INCY)
#
#       .. Scalar Arguments ..
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
# > CDOTC forms the dot product of two complex vectors
# >      CDOTC = X^H * Y
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
# > \param[in] CY
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
# >     jack dongarra, linpack,  3/11/78.
# >     modified 12/3/93, array(1) declarations changed to array(*)
# > \endverbatim
# >
#  =====================================================================
def cdotc(N, CX, INCX, CY, INCY):
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
    # COMPLEX CX(*),CY(*)
    #     ..
    #
    #  =====================================================================
    CTEMP = 0
    if N <= 0:
        return 0
    if INCX == 1 and INCY == 1:
        # code for both increments equal to 1
        for I in range(N):
            CTEMP += CX[I].conjugate() * CY[I]
    else:
        # code for unequal increments or equal increments not equal to 1
        IX = 1
        IY = 1
        if INCX < 0:
            IX = (-N + 1) * INCX + 1
        if INCY < 0:
            IY = (-N + 1) * INCY + 1
        for I in range(N):
            CTEMP += CX[IX].conjugate() * CY[IY]
            IX += INCX
            IY += INCY
    return CTEMP
