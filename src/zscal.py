# > \brief \b ZSCAL
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def ZSCAL(N,ZA,ZX,INCX)
#
#       .. Scalar Arguments ..
#       COMPLEX*16 ZA
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
# >    ZSCAL scales a vector by a constant.
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
# > \param[in,out] ZX
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
# > \ingroup complex16_blas_level1
#
# > \par Further Details:
#  =====================
# >
# > \verbatim
# >
# >     jack dongarra, 3/11/78.
# >     modified 3/93 to return if incx <= 0.
# >     modified 12/3/93, array(1) declarations changed to array(*)
# > \endverbatim
# >
#  =====================================================================
def zscal(N, ZA, ZX, INCX):
    #
    #  -- Reference BLAS level1 routine (version 3.8.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     November 2017
    #
    #     .. Scalar Arguments ..
    # COMPLEX*16 ZA
    # INTEGER INCX,N
    #     ..
    #     .. Array Arguments ..
    # COMPLEX*16 ZX(*)
    #     ..
    #
    #  =====================================================================
    if N <= 0 or INCX <= 0:
        return
    if INCX == 1:
        # code for increment equal to 1
        for I in range(N):
            ZX[I] = ZA * ZX[I]
    else:
        # code for increment not equal to 1
        for I in range(0, N * INCX, INCX):
            ZX[I] = ZA * ZX[I]
