# > \brief \b SASUM
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def SASUM(N,SX,INCX)
#
#       .. Scalar Arguments ..
#       INTEGER INCX,N
#       ..
#       .. Array Arguments ..
#       REAL SX(*)
#       ..
#
#
# > \par Purpose:
#  =============
# >
# > \verbatim
# >
# >    SASUM takes the sum of the absolute values.
# >    uses unrolled loops for increment equal to one.
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
# > \param[in] SX
# > \verbatim
# >          SX is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
# > \endverbatim
# >
# > \param[in] INCX
# > \verbatim
# >          INCX is INTEGER
# >         storage spacing between elements of SX
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
# >     modified 3/93 to return if incx <= 0.
# >     modified 12/3/93, array(1) declarations changed to array(*)
# > \endverbatim
# >
#  =====================================================================
def SASUM(N, SX, INCX):
    #
    #  -- Reference BLAS level1 routine (version 3.8.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     November 2017
    #
    #     .. Scalar Arguments ..
    #   INTEGER INCX,N
    #     ..
    #     .. Array Arguments ..
    #   REAL SX(*)
    #     ..
    #
    #  =====================================================================
    #
    #     .. Local Scalars ..
    #   REAL STEMP
    #   INTEGER I,M,MP1,NINCX
    #     ..
    #     .. Intrinsic Functions ..
    #   INTRINSIC ABS,MOD
    #     ..
    STEMP = 0
    if N <= 0 or INCX <= 0:
        return
    if INCX == 1:
        #        code for increment equal to 1
        #
        #
        #        clean-up loop
        #
        M = N % 6
        if M != 0:
            for I in range(M):
                STEMP += abs(SX[I])
            if N < 6:
                return STEMP
        for I in range(M, N, 6):
            STEMP += (
                abs(SX[I])
                + abs(SX(I + 1))
                + abs(SX(I + 2))
                + abs(SX(I + 3))
                + abs(SX(I + 4))
                + abs(SX(I + 5))
            )
    else:
        # code for increment not equal to 1
        for I in range(0, N * INCX, INCX):
            STEMP += abs(SX[I])
    return STEMP
