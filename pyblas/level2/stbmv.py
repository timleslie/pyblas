# > \brief \b STBMV
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def STBMV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
#
#       .. Scalar Arguments ..
#       INTEGER INCX,K,LDA,N
#       CHARACTER DIAG,TRANS,UPLO
#       ..
#       .. Array Arguments ..
#       REAL A(LDA,*),X(*)
#       ..
#
#
# > \par Purpose:
#  =============
# >
# > \verbatim
# >
# > STBMV  performs one of the matrix-vector operations
# >
# >    x := A*x,   or   x := A**T*x,
# >
# > where x is an n element vector and  A is an n by n unit, or non-unit,
# > upper or lower triangular band matrix, with ( k + 1 ) diagonals.
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
# >           On entry, TRANS specifies the operation to be performed as
# >           follows:
# >
# >              TRANS = 'N' or 'n'   x := A*x.
# >
# >              TRANS = 'T' or 't'   x := A**T*x.
# >
# >              TRANS = 'C' or 'c'   x := A**T*x.
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
# > \param[in] K
# > \verbatim
# >          K is INTEGER
# >           On entry with UPLO = 'U' or 'u', K specifies the number of
# >           super-diagonals of the matrix A.
# >           On entry with UPLO = 'L' or 'l', K specifies the number of
# >           sub-diagonals of the matrix A.
# >           K must satisfy  0 <= K.
# > \endverbatim
# >
# > \param[in] A
# > \verbatim
# >          A is REAL array, dimension ( LDA, N )
# >           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
# >           by n part of the array A must contain the upper triangular
# >           band part of the matrix of coefficients, supplied column by
# >           column, with the leading diagonal of the matrix in row
# >           ( k + 1 ) of the array, the first super-diagonal starting at
# >           position 2 in row k, and so on. The top left k by k triangle
# >           of the array A is not referenced.
# >           The following program segment will transfer an upper
# >           triangular band matrix from conventional full matrix storage
# >           to band storage:
# >
# >                 DO 20, J = 1, N
# >                    M = K + 1 - J
# >                    DO 10, I = max( 1, J - K ), J
# >                       A( M + I, J ) = matrix( I, J )
# >              10    CONTINUE
# >              20 CONTINUE
# >
# >           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
# >           by n part of the array A must contain the lower triangular
# >           band part of the matrix of coefficients, supplied column by
# >           column, with the leading diagonal of the matrix in row 1 of
# >           the array, the first sub-diagonal starting at position 1 in
# >           row 2, and so on. The bottom right k by k triangle of the
# >           array A is not referenced.
# >           The following program segment will transfer a lower
# >           triangular band matrix from conventional full matrix storage
# >           to band storage:
# >
# >                 DO 20, J = 1, N
# >                    M = 1 - J
# >                    DO 10, I = J, min( N, J + K )
# >                       A( M + I, J ) = matrix( I, J )
# >              10    CONTINUE
# >              20 CONTINUE
# >
# >           Note that when DIAG = 'U' or 'u' the elements of the array A
# >           corresponding to the diagonal elements of the matrix are not
# >           referenced, but are assumed to be unity.
# > \endverbatim
# >
# > \param[in] LDA
# > \verbatim
# >          LDA is INTEGER
# >           On entry, LDA specifies the first dimension of A as declared
# >           in the calling (sub) program. LDA must be at least
# >           ( k + 1 ).
# > \endverbatim
# >
# > \param[in,out] X
# > \verbatim
# >          X is REAL array, dimension at least
# >           ( 1 + ( n - 1 )*abs( INCX ) ).
# >           Before entry, the incremented array X must contain the n
# >           element vector x. On exit, X is overwritten with the
# >           transformed vector x.
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
# > \ingroup single_blas_level2
#
# > \par Further Details:
#  =====================
# >
# > \verbatim
# >
# >  Level 2 Blas routine.
# >  The vector and matrix arguments are not referenced when N = 0, or M = 0
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


def STBMV(UPLO, TRANS, DIAG, N, K, A, LDA, X, INCX):
    #
    #  -- Reference BLAS level2 routine (version 3.7.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     December 2016
    #
    #     .. Scalar Arguments ..
    #   INTEGER INCX,K,LDA,N
    #   CHARACTER DIAG,TRANS,UPLO
    #     ..
    #     .. Array Arguments ..
    #   REAL A(LDA,*),X(*)
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
    #   INTEGER I,INFO,IX,J,JX,KPLUS1,KX,L
    #   LOGICAL NOUNIT
    #     ..
    #     .. External Functions ..
    #   LOGICAL LSAME
    #   EXTERNAL LSAME
    #     ..
    #     .. External Subroutines ..
    #   EXTERNAL XERBLA
    #     ..
    #     .. Intrinsic Functions ..
    #   INTRINSIC MAX,MIN
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
    elif K < 0:
        INFO = 5
    elif LDA < (K + 1):
        INFO = 7
    elif INCX == 0:
        INFO = 9
    if INFO != 0:
        xerbla("STBMV ", INFO)

    # Quick return if possible.
    if N == 0:
        return
    #
    NOUNIT = lsame(DIAG, "N")
    #
    #     Set up the start point in X if the increment is not unity. This
    #     will be  ( N - 1 )*INCX   too small for descending loops.
    #
    if INCX <= 0:
        KX = 1 - (N - 1) * INCX
    elif INCX != 1:
        KX = 1
    #
    #     Start the operations. In this version the elements of A are
    #     accessed sequentially with one pass through A.
    #
    if lsame(TRANS, "N"):
        # Form  x := A*x.
        if lsame(UPLO, "U"):
            KPLUS1 = K + 1
            if INCX == 1:
                for J in range(N):
                    if X[J] != 0:
                        TEMP = X[J]
                        L = KPLUS1 - J
                        for I in range(max(1, J - K) - 1, J - 1):
                            X[I] += TEMP * A[L + I, J]
                        if NOUNIT:
                            X[J] = X[J] * A[KPLUS1, J]
            else:
                JX = KX
                for J in range(N):
                    if X[JX] != 0:
                        TEMP = X[JX]
                        IX = KX
                        L = KPLUS1 - J
                        for I in range(max(1, J - K) - 1, J - 1):
                            X[IX] = X[IX] + TEMP * A[L + I, J]
                            IX += INCX
                        if NOUNIT:
                            X[JX] = X[JX] * A[KPLUS1, J]
                    JX += INCX
                    if J > K:
                        KX += INCX
        else:
            if INCX == 1:
                for J in range(N - 1, -1, -1):
                    if X[J] != 0:
                        TEMP = X[J]
                        L = 1 - J
                        for I in range(min(N, J + K) - 1, J - 1, -1):
                            X[I] += TEMP * A[L + I, J]
                        if NOUNIT:
                            X[J] = X[J] * A[1, J]
            else:
                KX += (N - 1) * INCX
                JX = KX
                for J in range(N - 1, -1, -1):
                    if X[JX] != 0:
                        TEMP = X[JX]
                        IX = KX
                        L = 1 - J
                        for I in range(min(N, J + K) - 1, J - 1, -1):
                            X[IX] = X[IX] + TEMP * A[L + I, J]
                            IX -= INCX
                        if NOUNIT:
                            X[JX] = X[JX] * A[1, J]
                    JX -= INCX
                    if (N - J) >= K:
                        KX -= INCX
    else:
        #
        #        Form  x := A**T*x.
        #
        if lsame(UPLO, "U"):
            KPLUS1 = K + 1
            if INCX == 1:
                for J in range(N - 1, -1, -1):
                    TEMP = X[J]
                    L = KPLUS1 - J
                    if NOUNIT:
                        TEMP = TEMP * A[KPLUS1, J]
                    for I in range(J - 2, max(1, J - K) - 2, -1):
                        TEMP += A[L + I, J] * X[I]
                    X[J] = TEMP
            else:
                KX += (N - 1) * INCX
                JX = KX
                for J in range(N - 1, -1, -1):
                    TEMP = X[JX]
                    KX -= INCX
                    IX = KX
                    L = KPLUS1 - J
                    if NOUNIT:
                        TEMP = TEMP * A[KPLUS1, J]
                    for I in range(J - 2, max(1, J - K) - 2, -1):
                        TEMP += A[L + I, J] * X[IX]
                        IX -= INCX
                    X[JX] = TEMP
                    JX -= INCX
        else:
            if INCX == 1:
                for J in range(N):
                    TEMP = X[J]
                    L = 1 - J
                    if NOUNIT:
                        TEMP = TEMP * A[1, J]
                    for I in range(J, min(N, J + K)):
                        TEMP += A[L + I, J] * X[I]
                    X[J] = TEMP
            else:
                JX = KX
                for J in range(N):
                    TEMP = X[JX]
                    KX += INCX
                    IX = KX
                    L = 1 - J
                    if NOUNIT:
                        TEMP = TEMP * A[1, J]
                    for I in range(J, min(N, J + K)):
                        TEMP += A[L + I, J] * X[IX]
                        IX += INCX
                    X[JX] = TEMP
                    JX += INCX
