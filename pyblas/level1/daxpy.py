# > \brief \b DAXPY
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def DAXPY(N,DA,DX,INCX,DY,INCY)
#
#       .. Scalar Arguments ..
#       DOUBLE PRECISION DA
#       INTEGER INCX,INCY,N
#       ..
#       .. Array Arguments ..
#       DOUBLE PRECISION DX(*),DY(*)
#       ..
#
#
# > \par Purpose:
#  =============
# >
# > \verbatim
# >
# >    DAXPY constant times a vector plus a vector.
# >    uses unrolled loops for increments equal to one.
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
# >
# > \param[in,out] DY
# > \verbatim
# >          DY is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
# > \endverbatim
# >
# > \param[in] INCY
# > \verbatim
# >          INCY is INTEGER
# >         storage spacing between elements of DY
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
# >     modified 12/3/93, array(1) declarations changed to array(*)
# > \endverbatim
# >
#  =====================================================================
def DAXPY(N, DA, DX, INCX, DY, INCY):
    #
    #  -- Reference BLAS level1 routine (version 3.8.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     November 2017
    #
    #     .. Scalar Arguments ..
    # DOUBLE PRECISION DA
    # INTEGER INCX,INCY,N
    #     ..
    #     .. Array Arguments ..
    # DOUBLE PRECISION DX(*),DY(*)
    #     ..
    #
    #  =====================================================================
    #
    #     .. Local Scalars ..
    # INTEGER I,IX,IY,M,MP1
    #     ..
    #     .. Intrinsic Functions ..
    # INTRINSIC MOD
    #     ..
    if N <= 0:
        return
    if DA == 0.0:
        return
    if INCX == 1 and INCY == 1:
        #
        # code for both increments equal to 1
        #
        #
        # clean-up loop
        #
        M = N % 4
        if M != 0:
            for I in range(M):
                DY[I] += DA * DX[I]
        if N < 4:
            return
        for I in range(M, N, 4):
            DY[I] += DA * DX[I]
            DY[I + 1] += DA * DX[I + 1]
            DY[I + 2] += DA * DX[I + 2]
            DY[I + 3] += DA * DX[I + 3]
    else:
        # code for unequal increments or equal increments not equal to 1
        IX = 1
        IY = 1
        if INCX < 0:
            IX = (-N + 1) * INCX + 1
        if INCY < 0:
            IY = (-N + 1) * INCY + 1
        for I in range(N):
            DY[IY] += DA * DX[IX]
            IX += INCX
            IY += INCY
    return
