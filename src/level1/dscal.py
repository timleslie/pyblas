# > \brief \b DSCAL
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def DSCAL(N,DA,DX,INCX)
#
#       .. Scalar Arguments ..
#       DOUBLE PRECISION DA
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
# >    DSCAL scales a vector by a constant.
# >    uses unrolled loops for increment equal to 1.
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
# > \param[in] DA
# > \verbatim
# >          DA is DOUBLE PRECISION
# >           On entry, DA specifies the scalar alpha.
# > \endverbatim
# >
# > \param[in,out] DX
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
def DSCAL(N, DA, DX, INCX):
    #
    #  -- Reference BLAS level1 routine (version 3.8.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     November 2017
    #
    #     .. Scalar Arguments ..
    # DOUBLE PRECISION DA
    # INTEGER INCX,N
    #     ..
    #     .. Array Arguments ..
    # DOUBLE PRECISION DX(*)
    #     ..
    #
    #  =====================================================================
    #
    #     .. Local Scalars ..
    # INTEGER I,M,MP1,NINCX
    #     ..
    #     .. Intrinsic Functions ..
    # INTRINSIC MOD
    #     ..
    if N <= 0 or INCX <= 0:
        return
    if INCX == 1:
        #
        #        code for increment equal to 1
        #
        #
        #        clean-up loop
        #
        M = N % 5
        if M != 0:
            for I in range(M):
                DX[I] = DA * DX[I]
            if N < 5:
                return
        for I in range(M, N, 5):
            DX[I] = DA * DX[I]
            DX[I + 1] = DA * DX[I + 1]
            DX[I + 2] = DA * DX[I + 2]
            DX[I + 3] = DA * DX[I + 3]
            DX[I + 4] = DA * DX[I + 4]
    else:
        #
        #        code for increment not equal to 1
        #
        for I in range(0, N * INCX, INCX):
            DX[I] = DA * DX[I]
