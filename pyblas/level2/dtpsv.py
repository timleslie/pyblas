# > \brief \b DTPSV
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def DTPSV(UPLO,TRANS,DIAG,N,AP,X,INCX)
#
#       .. Scalar Arguments ..
#       INTEGER INCX,N
#       CHARACTER DIAG,TRANS,UPLO
#       ..
#       .. Array Arguments ..
#       DOUBLE PRECISION AP(*),X(*)
#       ..
#
#
# > \par Purpose:
#  =============
# >
# > \verbatim
# >
# > DTPSV  solves one of the systems of equations
# >
# >    A*x = b,   or   A**T*x = b,
# >
# > where b and x are n element vectors and A is an n by n unit, or
# > non-unit, upper or lower triangular matrix, supplied in packed form.
# >
# > No test for singularity or near-singularity is included in this
# > routine. Such tests must be performed before calling this routine.
# > \endverbatim
#
#  Arguments:
#  ==========
#
# > \param[in] UPLO
# > \verbatim
# >          UPLO is CHARACTER*1
# >           On entry, UPLO specifies whether the matrix is an upper or
# >           lower triangular matrix as follows:
# >
# >              UPLO = 'U' or 'u'   A is an upper triangular matrix.
# >
# >              UPLO = 'L' or 'l'   A is a lower triangular matrix.
# > \endverbatim
# >
# > \param[in] TRANS
# > \verbatim
# >          TRANS is CHARACTER*1
# >           On entry, TRANS specifies the equations to be solved as
# >           follows:
# >
# >              TRANS = 'N' or 'n'   A*x = b.
# >
# >              TRANS = 'T' or 't'   A**T*x = b.
# >
# >              TRANS = 'C' or 'c'   A**T*x = b.
# > \endverbatim
# >
# > \param[in] DIAG
# > \verbatim
# >          DIAG is CHARACTER*1
# >           On entry, DIAG specifies whether or not A is unit
# >           triangular as follows:
# >
# >              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
# >
# >              DIAG = 'N' or 'n'   A is not assumed to be unit
# >                                  triangular.
# > \endverbatim
# >
# > \param[in] N
# > \verbatim
# >          N is INTEGER
# >           On entry, N specifies the order of the matrix A.
# >           N must be at least zero.
# > \endverbatim
# >
# > \param[in] AP
# > \verbatim
# >          AP is DOUBLE PRECISION array, dimension at least
# >           ( ( n*( n + 1 ) )/2 ).
# >           Before entry with  UPLO = 'U' or 'u', the array AP must
# >           contain the upper triangular matrix packed sequentially,
# >           column by column, so that AP( 1 ) contains a( 1, 1 ),
# >           AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 )
# >           respectively, and so on.
# >           Before entry with UPLO = 'L' or 'l', the array AP must
# >           contain the lower triangular matrix packed sequentially,
# >           column by column, so that AP( 1 ) contains a( 1, 1 ),
# >           AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 )
# >           respectively, and so on.
# >           Note that when  DIAG = 'U' or 'u', the diagonal elements of
# >           A are not referenced, but are assumed to be unity.
# > \endverbatim
# >
# > \param[in,out] X
# > \verbatim
# >          X is DOUBLE PRECISION array, dimension at least
# >           ( 1 + ( n - 1 )*abs( INCX ) ).
# >           Before entry, the incremented array X must contain the n
# >           element right-hand side vector b. On exit, X is overwritten
# >           with the solution vector x.
# > \endverbatim
# >
# > \param[in] INCX
# > \verbatim
# >          INCX is INTEGER
# >           On entry, INCX specifies the increment for the elements of
# >           X. INCX must not be zero.
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
# > \ingroup double_blas_level2
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


def DTPSV(UPLO, TRANS, DIAG, N, AP, X, INCX):
    #
    #  -- Reference BLAS level2 routine (version 3.7.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     December 2016
    #
    #     .. Scalar Arguments ..
    #   INTEGER INCX,N
    #   CHARACTER DIAG,TRANS,UPLO
    #     ..
    #     .. Array Arguments ..
    #   DOUBLE PRECISION AP(*),X(*)
    #     ..
    #
    #  =====================================================================
    #
    #     .. Parameters ..
    #   DOUBLE PRECISION ZERO
    #   PARAMETER (ZERO=0.0D+0)
    #     ..
    #     .. Local Scalars ..
    #   DOUBLE PRECISION TEMP
    #   INTEGER I,INFO,IX,J,JX,K,KK,KX
    #   LOGICAL NOUNIT
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
    elif not lsame(TRANS, "N") and not lsame(TRANS, "T") and not lsame(TRANS, "C"):
        INFO = 2
    elif not lsame(DIAG, "U") and not lsame(DIAG, "N"):
        INFO = 3
    elif N < 0:
        INFO = 4
    elif INCX == 0:
        INFO = 7
    if INFO != 0:
        xerbla("DTPSV ", INFO)

    # Quick return if possible.
    if N == 0:
        return
    #
    NOUNIT = lsame(DIAG, "N")
    #
    #     Set up the start point in X if the increment is not unity. This
    #     will be  ( N - 1 )*INCX  too small for descending loops.
    #
    if INCX <= 0:
        KX = 1 - (N - 1) * INCX
    elif INCX != 1:
        KX = 1
    #
    #     Start the operations. In this version the elements of AP are
    #     accessed sequentially with one pass through AP.
    #
    if lsame(TRANS, "N"):
        #
        #        Form  x := inv( A )*x.
        #
        if lsame(UPLO, "U"):
            KK = (N * (N + 1)) / 2
            if INCX == 1:
                for J in range(N - 1, -1, -1):
                    if X[J] != 0:
                        if NOUNIT:
                            X[J] = X[J] / AP[KK]
                        TEMP = X[J]
                        K = KK - 1
                        for I in range(J - 2, -1, -1):
                            X[I] -= TEMP * AP[K]
                            K -= 1
                    KK -= J
            else:
                JX = KX + (N - 1) * INCX
                for J in range(N - 1, -1, -1):
                    if X[JX] != 0:
                        if NOUNIT:
                            X[JX] = X[JX] / AP[KK]
                        TEMP = X[JX]
                        IX = JX
                        for K in range(KK - 2, KK - J - 2, -1):
                            IX -= INCX
                            X[IX] -= TEMP * AP[K]
                    JX -= INCX
                    KK -= J
        else:
            KK = 1
            if INCX == 1:
                for J in range(N):
                    if X[J] != 0:
                        if NOUNIT:
                            X[J] = X[J] / AP[KK]
                        TEMP = X[J]
                        K = KK + 1
                        for I in range(J, N):
                            X[I] -= TEMP * AP[K]
                            K += 1
                    KK += N - J + 1
            else:
                JX = KX
                for J in range(N):
                    if X[JX] != 0:
                        if NOUNIT:
                            X[JX] = X[JX] / AP[KK]
                        TEMP = X[JX]
                        IX = JX
                        for K in range(KK, KK + N - J):
                            IX += INCX
                            X[IX] -= TEMP * AP[K]
                    JX += INCX
                    KK += N - J + 1
    else:
        #
        #        Form  x := inv( A**T )*x.
        #
        if lsame(UPLO, "U"):
            KK = 1
            if INCX == 1:
                for J in range(N):
                    TEMP = X[J]
                    K = KK
                    for I in range(J - 1):
                        TEMP -= AP[K] * X[I]
                        K += 1
                    if NOUNIT:
                        TEMP = TEMP / AP[KK + J - 1]
                    X[J] = TEMP
                    KK += J
            else:
                JX = KX
                for J in range(N):
                    TEMP = X[JX]
                    IX = KX
                    for K in range(KK - 1, KK + J - 2):
                        TEMP -= AP[K] * X[IX]
                        IX += INCX
                    if NOUNIT:
                        TEMP = TEMP / AP[KK + J - 1]
                    X[JX] = TEMP
                    JX += INCX
                    KK += J
        else:
            KK = (N * (N + 1)) / 2
            if INCX == 1:
                for J in range(N - 1, -1, -1):
                    TEMP = X[J]
                    K = KK
                    for I in range(N - 1, J - 1, -1):
                        TEMP -= AP[K] * X[I]
                        K -= 1
                    if NOUNIT:
                        TEMP = TEMP / AP[KK - N + J]
                    X[J] = TEMP
                    KK -= N - J + 1
            else:
                KX += (N - 1) * INCX
                JX = KX
                for J in range(N - 1, -1, -1):
                    TEMP = X[JX]
                    IX = KX
                    for K in range(KK - 1, KK - (N - (J + 1)) - 2, -1):
                        TEMP -= AP[K] * X[IX]
                        IX -= INCX
                    if NOUNIT:
                        TEMP /= AP[KK - N + J]
                    X[JX] = TEMP
                    JX -= INCX
                    KK -= N - J + 1
