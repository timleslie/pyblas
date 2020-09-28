# > \brief \b IDAMAX
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       INTEGER FUNCTION IDAmax(N,DX,INCX)
#
#       .. Scalar Arguments ..
#       INTEGER INCX,N
#       ..
#       .. Array Arguments ..
#       DOUBLE PRECISION DX(*)
#       ..
#
#
# > \par Purpose:
#  =============
# >
# > \verbatim
# >
# >    IDAMAX finds the index of the first element having maximum absolute value.
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
# > \param[in] DX
# > \verbatim
# >          DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
# > \endverbatim
# >
# > \param[in] INCX
# > \verbatim
# >          INCX is INTEGER
# >         storage spacing between elements of DX
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
def IDAmax(N, DX, INCX):
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
    # DOUBLE PRECISION DX(*)
    #     ..
    #
    #  =====================================================================

    if N < 1 or INCX <= 0:
        return 0
    if N == 1:
        return 1

    IDAMAX = 1
    if INCX == 1:
        #
        #        code for increment equal to 1
        #
        DMAX = abs(DX(1))
        for I in range(1, N):
            if abs(DX[I]) > DMAX:
                IDAMAX = I
                DMAX = abs(DX[I])
    else:
        #
        #        code for increment not equal to 1
        #
        IX = 1
        DMAX = abs(DX(1))
        IX += INCX
        for I in range(1, N):
            if abs(DX(IX)) > DMAX:
                IDAMAX = I
                DMAX = abs(DX(IX))
            IX += INCX
    return IDAMAX
