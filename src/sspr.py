# > \brief \b SSPR
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def SSPR(UPLO,N,ALPHA,X,INCX,AP)
#
#       .. Scalar Arguments ..
#       REAL ALPHA
#       INTEGER INCX,N
#       CHARACTER UPLO
#       ..
#       .. Array Arguments ..
#       REAL AP(*),X(*)
#       ..
#
#
# > \par Purpose:
#  =============
# >
# > \verbatim
# >
# > SSPR    performs the symmetric rank 1 operation
# >
# >    A := alpha*x*x**T + A,
# >
# > where alpha is a real scalar, x is an n element vector and A is an
# > n by n symmetric matrix, supplied in packed form.
# > \endverbatim
#
#  Arguments:
#  ==========
#
# > \param[in] UPLO
# > \verbatim
# >          UPLO is CHARACTER*1
# >           On entry, UPLO specifies whether the upper or lower
# >           triangular part of the matrix A is supplied in the packed
# >           array AP as follows:
# >
# >              UPLO = 'U' or 'u'   The upper triangular part of A is
# >                                  supplied in AP.
# >
# >              UPLO = 'L' or 'l'   The lower triangular part of A is
# >                                  supplied in AP.
# > \endverbatim
# >
# > \param[in] N
# > \verbatim
# >          N is INTEGER
# >           On entry, N specifies the order of the matrix A.
# >           N must be at least zero.
# > \endverbatim
# >
# > \param[in] ALPHA
# > \verbatim
# >          ALPHA is REAL
# >           On entry, ALPHA specifies the scalar alpha.
# > \endverbatim
# >
# > \param[in] X
# > \verbatim
# >          X is REAL array, dimension at least
# >           ( 1 + ( n - 1 )*abs( INCX ) ).
# >           Before entry, the incremented array X must contain the n
# >           element vector x.
# > \endverbatim
# >
# > \param[in] INCX
# > \verbatim
# >          INCX is INTEGER
# >           On entry, INCX specifies the increment for the elements of
# >           X. INCX must not be zero.
# > \endverbatim
# >
# > \param[in,out] AP
# > \verbatim
# >          AP is REAL array, dimension at least
# >           ( ( n*( n + 1 ) )/2 ).
# >           Before entry with  UPLO = 'U' or 'u', the array AP must
# >           contain the upper triangular part of the symmetric matrix
# >           packed sequentially, column by column, so that AP( 1 )
# >           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
# >           and a( 2, 2 ) respectively, and so on. On exit, the array
# >           AP is overwritten by the upper triangular part of the
# >           updated matrix.
# >           Before entry with UPLO = 'L' or 'l', the array AP must
# >           contain the lower triangular part of the symmetric matrix
# >           packed sequentially, column by column, so that AP( 1 )
# >           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
# >           and a( 3, 1 ) respectively, and so on. On exit, the array
# >           AP is overwritten by the lower triangular part of the
# >           updated matrix.
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
# > \date December 2016
#
# > \ingroup single_blas_level2
#
# > \par Further Details:
#  =====================
# >
# > \verbatim
# >
# >  Level 2 Blas routine.
# >
# >  -- Written on 22-October-1986.
# >     Jack Dongarra, Argonne National Lab.
# >     Jeremy Du Croz, Nag Central Office.
# >     Sven Hammarling, Nag Central Office.
# >     Richard Hanson, Sandia National Labs.
# > \endverbatim
# >
#  =====================================================================
from util import lsame
from xerbla import xerbla


def SSPR(UPLO, N, ALPHA, X, INCX, AP):
    #
    #  -- Reference BLAS level2 routine (version 3.7.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     December 2016
    #
    #     .. Scalar Arguments ..
    #   REAL ALPHA
    #   INTEGER INCX,N
    #   CHARACTER UPLO
    #     ..
    #     .. Array Arguments ..
    #   REAL AP(*),X(*)
    #     ..
    #
    #  =====================================================================
    #
    #     .. Parameters ..
    #   REAL ZERO
    #   PARAMETER (ZERO=0.0E+0)
    #     ..
    #     .. Local Scalars ..
    #   REAL TEMP
    #   INTEGER I,INFO,IX,J,JX,K,KK,KX
    #     ..
    #     .. External Functions ..
    #   LOGICAL LSAME
    #   EXTERNAL LSAME
    #     ..
    #     .. External Subroutines ..
    #   EXTERNAL XERBLA
    #     ..

    # Test the input parameters.
    INFO = 0
    if not lsame(UPLO, "U") and not lsame(UPLO, "L"):
        INFO = 1
    elif N < 0:
        INFO = 2
    elif INCX == 0:
        INFO = 5
    if INFO != 0:
        xerbla("SSPR  ", INFO)

    # Quick return if possible.
    if (N == 0) or (ALPHA == 0):
        return
    #
    #     Set the start point in X if the increment is not unity.
    #
    if INCX <= 0:
        KX = 1 - (N - 1) * INCX
    elif INCX != 1:
        KX = 1
    #
    #     Start the operations. In this version the elements of the array AP
    #     are accessed sequentially with one pass through AP.
    #
    KK = 1
    if lsame(UPLO, "U"):
        #
        #        Form  A  when upper triangle is stored in AP.
        #
        if INCX == 1:
            for J in range(N):
                if X[J] != 0:
                    TEMP = ALPHA * X[J]
                    K = KK
                    for I in range(J):
                        AP[K] = AP[K] + X[I] * TEMP
                        K += 1
                KK += J
        else:
            JX = KX
            for J in range(N):
                if X[JX] != 0:
                    TEMP = ALPHA * X[JX]
                    IX = KX
                    for K in range(KK - 1, KK + J - 1):
                        AP[K] = AP[K] + X[IX] * TEMP
                        IX += INCX
                JX += INCX
                KK += J
    else:
        #
        #        Form  A  when lower triangle is stored in AP.
        #
        if INCX == 1:
            for J in range(N):
                if X[J] != 0:
                    TEMP = ALPHA * X[J]
                    K = KK
                    for I in range(J - 1, N):
                        AP[K] = AP[K] + X[I] * TEMP
                        K += 1
                KK = KK + N - J + 1
        else:
            JX = KX
            for J in range(N):
                if X[JX] != 0:
                    TEMP = ALPHA * X[JX]
                    IX = JX
                    for K in range(KK - 1, KK + N - J):
                        AP[K] = AP[K] + X[IX] * TEMP
                        IX += INCX
                JX += INCX
                KK = KK + N - J + 1
