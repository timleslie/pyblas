# > \brief \b CHBMV
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def CHBMV(UPLO,N,K,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
#
#       .. Scalar Arguments ..
#       COMPLEX ALPHA,BETA
#       INTEGER INCX,INCY,K,LDA,N
#       CHARACTER UPLO
#       ..
#       .. Array Arguments ..
#       COMPLEX A(LDA,*),X(*),Y(*)
#       ..
#
#
# > \par Purpose:
#  =============
# >
# > \verbatim
# >
# > CHBMV  performs the matrix-vector  operation
# >
# >    y := alpha*A*x + beta*y,
# >
# > where alpha and beta are scalars, x and y are n element vectors and
# > A is an n by n hermitian band matrix, with k super-diagonals.
# > \endverbatim
#
#  Arguments:
#  ==========
#
# > \param[in] UPLO
# > \verbatim
# >          UPLO is CHARACTER*1
# >           On entry, UPLO specifies whether the upper or lower
# >           triangular part of the band matrix A is being supplied as
# >           follows:
# >
# >              UPLO = 'U' or 'u'   The upper triangular part of A is
# >                                  being supplied.
# >
# >              UPLO = 'L' or 'l'   The lower triangular part of A is
# >                                  being supplied.
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
# >           On entry, K specifies the number of super-diagonals of the
# >           matrix A. K must satisfy  0 <= K.
# > \endverbatim
# >
# > \param[in] ALPHA
# > \verbatim
# >          ALPHA is COMPLEX
# >           On entry, ALPHA specifies the scalar alpha.
# > \endverbatim
# >
# > \param[in] A
# > \verbatim
# >          A is COMPLEX array, dimension ( LDA, N )
# >           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
# >           by n part of the array A must contain the upper triangular
# >           band part of the hermitian matrix, supplied column by
# >           column, with the leading diagonal of the matrix in row
# >           ( k + 1 ) of the array, the first super-diagonal starting at
# >           position 2 in row k, and so on. The top left k by k triangle
# >           of the array A is not referenced.
# >           The following program segment will transfer the upper
# >           triangular part of a hermitian band matrix from conventional
# >           full matrix storage to band storage:
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
# >           band part of the hermitian matrix, supplied column by
# >           column, with the leading diagonal of the matrix in row 1 of
# >           the array, the first sub-diagonal starting at position 1 in
# >           row 2, and so on. The bottom right k by k triangle of the
# >           array A is not referenced.
# >           The following program segment will transfer the lower
# >           triangular part of a hermitian band matrix from conventional
# >           full matrix storage to band storage:
# >
# >                 DO 20, J = 1, N
# >                    M = 1 - J
# >                    DO 10, I = J, min( N, J + K )
# >                       A( M + I, J ) = matrix( I, J )
# >              10    CONTINUE
# >              20 CONTINUE
# >
# >           Note that the imaginary parts of the diagonal elements need
# >           not be set and are assumed to be zero.
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
# > \param[in] X
# > \verbatim
# >          X is COMPLEX array, dimension at least
# >           ( 1 + ( n - 1 )*abs( INCX ) ).
# >           Before entry, the incremented array X must contain the
# >           vector x.
# > \endverbatim
# >
# > \param[in] INCX
# > \verbatim
# >          INCX is INTEGER
# >           On entry, INCX specifies the increment for the elements of
# >           X. INCX must not be zero.
# > \endverbatim
# >
# > \param[in] BETA
# > \verbatim
# >          BETA is COMPLEX
# >           On entry, BETA specifies the scalar beta.
# > \endverbatim
# >
# > \param[in,out] Y
# > \verbatim
# >          Y is COMPLEX array, dimension at least
# >           ( 1 + ( n - 1 )*abs( INCY ) ).
# >           Before entry, the incremented array Y must contain the
# >           vector y. On exit, Y is overwritten by the updated vector y.
# > \endverbatim
# >
# > \param[in] INCY
# > \verbatim
# >          INCY is INTEGER
# >           On entry, INCY specifies the increment for the elements of
# >           Y. INCY must not be zero.
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
# > \ingroup complex_blas_level2
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


def CHBMV(UPLO, N, K, ALPHA, A, LDA, X, INCX, BETA, Y, INCY):
    #
    #  -- Reference BLAS level2 routine (version 3.7.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     December 2016
    #
    #     .. Scalar Arguments ..
    #   COMPLEX ALPHA,BETA
    #   INTEGER INCX,INCY,K,LDA,N
    #   CHARACTER UPLO
    #     ..
    #     .. Array Arguments ..
    #   COMPLEX A(LDA,*),X(*),Y(*)
    #     ..
    #
    #  =====================================================================
    #
    #     .. Parameters ..
    #   COMPLEX ONE
    #   PARAMETER (ONE= (1.0E+0,0.0E+0))
    #   COMPLEX ZERO
    #   PARAMETER (ZERO= (0.0E+0,0.0E+0))
    #     ..
    #     .. Local Scalars ..
    #   COMPLEX TEMP1,TEMP2
    #   INTEGER I,INFO,IX,IY,J,JX,JY,KPLUS1,KX,KY,L
    #     ..
    #     .. External Functions ..
    #   LOGICAL LSAME
    #   EXTERNAL LSAME
    #     ..
    #     .. External Subroutines ..
    #   EXTERNAL XERBLA
    #     ..
    #     .. Intrinsic Functions ..
    #   INTRINSIC CONJG,MAX,MIN,REAL
    #     ..

    # Test the input parameters.
    INFO = 0
    if not lsame(UPLO, "U") and not lsame(UPLO, "L"):
        INFO = 1
    elif N < 0:
        INFO = 2
    elif K < 0:
        INFO = 3
    elif LDA < (K + 1):
        INFO = 6
    elif INCX == 0:
        INFO = 8
    elif INCY == 0:
        INFO = 11
    if INFO != 0:
        xerbla("CHBMV ", INFO)
        return

    # Quick return if possible.
    if (N == 0) or ((ALPHA == 0) and (BETA == 1)):
        return
    #
    #     Set up the start points in  X  and  Y.
    #
    if INCX > 0:
        KX = 1
    else:
        KX = 1 - (N - 1) * INCX
    if INCY > 0:
        KY = 1
    else:
        KY = 1 - (N - 1) * INCY
    #
    #     Start the operations. In this version the elements of the array A
    #     are accessed sequentially with one pass through A.
    #
    #     First form  y := beta*y.
    #
    if BETA != 1:
        if INCY == 1:
            if BETA == 0:
                for I in range(N):
                    Y[I] = 0
            else:
                for I in range(N):
                    Y[I] = BETA * Y[I]
        else:
            IY = KY
            if BETA == 0:
                for I in range(N):
                    Y[IY] = 0
                    IY += INCY
            else:
                for I in range(N):
                    Y[IY] = BETA * Y[IY]
                    IY += INCY
    if ALPHA == 0:
        return
    if lsame(UPLO, "U"):
        #
        #        Form  y  when upper triangle of A is stored.
        #
        KPLUS1 = K + 1
        if (INCX == 1) and (INCY == 1):
            for J in range(N):
                TEMP1 = ALPHA * X[J]
                TEMP2 = 0
                L = KPLUS1 - J
                for I in range(max(1, J - K) - 1, J - 1):
                    Y[I] += TEMP1 * A[L + I, J]
                    TEMP2 += A[L + I, J].conjugate() * X[I]
                Y[J] += TEMP1 * A[KPLUS1, J].real + ALPHA * TEMP2
        else:
            JX = KX
            JY = KY
            for J in range(N):
                TEMP1 = ALPHA * X[JX]
                TEMP2 = 0
                IX = KX
                IY = KY
                L = KPLUS1 - J
                for I in range(max(1, J - K) - 1, J - 1):
                    Y[IY] += TEMP1 * A[L + I, J]
                    TEMP2 += A[L + I, J].conjugate() * X[IX]
                    IX += INCX
                    IY += INCY
                Y[JY] += TEMP1 * A[KPLUS1, J].real + ALPHA * TEMP2
                JX += INCX
                JY += INCY
                if J > K:
                    KX += INCX
                    KY += INCY
    else:
        # Form  y  when lower triangle of A is stored.
        if (INCX == 1) and (INCY == 1):
            for J in range(N):
                TEMP1 = ALPHA * X[J]
                TEMP2 = 0
                Y[J] += TEMP1 * A[1, J].real
                L = 1 - J
                for I in range(J, min(N, J + K)):
                    Y[I] += TEMP1 * A[L + I, J]
                    TEMP2 += A[L + I, J].conjugate() * X[I]
                Y[J] += ALPHA * TEMP2
        else:
            JX = KX
            JY = KY
            for J in range(N):
                TEMP1 = ALPHA * X[JX]
                TEMP2 = 0
                Y[JY] += TEMP1 * A[1, J].real
                L = 1 - J
                IX = JX
                IY = JY
                for I in range(J, min(N, J + K)):
                    IX += INCX
                    IY += INCY
                    Y[IY] += TEMP1 * A[L + I, J]
                    TEMP2 += A[L + I, J].conjugate() * X[IX]
                Y[JY] += ALPHA * TEMP2
                JX += INCX
                JY += INCY
