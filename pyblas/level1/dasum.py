# > \brief \b DASUM
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def DASUM(N,DX,INCX)
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
# >    DASUM takes the sum of the absolute values.
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
# > \ingroup double_blas_level1
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
def dasum(N, DX, INCX):
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
    if N <= 0 or INCX <= 0:
        return 0
    DTEMP = 0
    if INCX == 1:
        # code for increment equal to 1

        # clean-up loop
        M = N % 6
        if M != 0:
            for I in range(M):
                DTEMP += abs(DX[I])
            if N < 6:
                return DTEMP
        for I in range(M, N, 6):
            DTEMP += (
                abs(DX[I])
                + abs(DX[I + 1])
                + abs(DX[I + 2])
                + abs(DX[I + 3])
                + abs(DX[I + 4])
                + abs(DX[I + 5])
            )
    else:
        # code for increment not equal to 1
        for I in range(0, N * INCX, INCX):
            DTEMP += abs(DX[I])
    return DTEMP
