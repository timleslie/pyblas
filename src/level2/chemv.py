# > \brief \b CHEMV
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def CHEMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
#
#       .. Scalar Arguments ..
#       COMPLEX ALPHA,BETA
#       INTEGER INCX,INCY,LDA,N
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
# > CHEMV  performs the matrix-vector  operation
# >
# >    y := alpha*A*x + beta*y,
# >
# > where alpha and beta are scalars, x and y are n element vectors and
# > A is an n by n hermitian matrix.
# > \endverbatim
#
#  Arguments:
#  ==========
#
# > \param[in] UPLO
# > \verbatim
# >          UPLO is CHARACTER*1
# >           On entry, UPLO specifies whether the upper or lower
# >           triangular part of the array A is to be referenced as
# >           follows:
# >
# >              UPLO = 'U' or 'u'   Only the upper triangular part of A
# >                                  is to be referenced.
# >
# >              UPLO = 'L' or 'l'   Only the lower triangular part of A
# >                                  is to be referenced.
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
# >          ALPHA is COMPLEX
# >           On entry, ALPHA specifies the scalar alpha.
# > \endverbatim
# >
# > \param[in] A
# > \verbatim
# >          A is COMPLEX array, dimension ( LDA, N )
# >           Before entry with  UPLO = 'U' or 'u', the leading n by n
# >           upper triangular part of the array A must contain the upper
# >           triangular part of the hermitian matrix and the strictly
# >           lower triangular part of A is not referenced.
# >           Before entry with UPLO = 'L' or 'l', the leading n by n
# >           lower triangular part of the array A must contain the lower
# >           triangular part of the hermitian matrix and the strictly
# >           upper triangular part of A is not referenced.
# >           Note that the imaginary parts of the diagonal elements need
# >           not be set and are assumed to be zero.
# > \endverbatim
# >
# > \param[in] LDA
# > \verbatim
# >          LDA is INTEGER
# >           On entry, LDA specifies the first dimension of A as declared
# >           in the calling (sub) program. LDA must be at least
# >           max( 1, n ).
# > \endverbatim
# >
# > \param[in] X
# > \verbatim
# >          X is COMPLEX array, dimension at least
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
# > \param[in] BETA
# > \verbatim
# >          BETA is COMPLEX
# >           On entry, BETA specifies the scalar beta. When BETA is
# >           supplied as zero then Y need not be set on input.
# > \endverbatim
# >
# > \param[in,out] Y
# > \verbatim
# >          Y is COMPLEX array, dimension at least
# >           ( 1 + ( n - 1 )*abs( INCY ) ).
# >           Before entry, the incremented array Y must contain the n
# >           element vector y. On exit, Y is overwritten by the updated
# >           vector y.
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


def CHEMV(UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY):
    #
    #  -- Reference BLAS level2 routine (version 3.7.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     December 2016
    #
    #     .. Scalar Arguments ..
    #   COMPLEX ALPHA,BETA
    #   INTEGER INCX,INCY,LDA,N
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
    #   INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY
    #     ..
    #     .. External Functions ..
    #   LOGICAL LSAME
    #   EXTERNAL LSAME
    #     ..
    #     .. External Subroutines ..
    #   EXTERNAL XERBLA
    #     ..
    #     .. Intrinsic Functions ..
    #   INTRINSIC CONJG,MAX,REAL
    #     ..

    # Test the input parameters.
    INFO = 0
    if not lsame(UPLO, "U") and not lsame(UPLO, "L"):
        INFO = 1
    elif N < 0:
        INFO = 2
    elif LDA < max(1, N):
        INFO = 5
    elif INCX == 0:
        INFO = 7
    elif INCY == 0:
        INFO = 10
    if INFO != 0:
        xerbla("CHEMV ", INFO)
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
    #     Start the operations. In this version the elements of A are
    #     accessed sequentially with one pass through the triangular part
    #     of A.
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
        #        Form  y  when A is stored in upper triangle.
        #
        if (INCX == 1) and (INCY == 1):
            for J in range(N):
                TEMP1 = ALPHA * X[J]
                TEMP2 = 0
                for I in range(J - 1):
                    Y[I] += TEMP1 * A[I, J]
                    TEMP2 += A[1, J].conjugate() * X[I]
                Y[J] += TEMP1 * A[J, J].real + ALPHA * TEMP2
        else:
            JX = KX
            JY = KY
            for J in range(N):
                TEMP1 = ALPHA * X[JX]
                TEMP2 = 0
                IX = KX
                IY = KY
                for I in range(J - 1):
                    Y[IY] += TEMP1 * A[I, J]
                    TEMP2 += A[1, J].conjugate() * X[IX]
                    IX += INCX
                    IY += INCY
                Y[JY] += TEMP1 * A[J, J].real + ALPHA * TEMP2
                JX += INCX
                JY += INCY
    else:
        #
        #        Form  y  when A is stored in lower triangle.
        #
        if (INCX == 1) and (INCY == 1):
            for J in range(N):
                TEMP1 = ALPHA * X[J]
                TEMP2 = 0
                Y[J] += TEMP1 * A[J, J].real
                for I in range(J, N):
                    Y[I] += TEMP1 * A[I, J]
                    TEMP2 += A[1, J].conjugate() * X[I]
                Y[J] += ALPHA * TEMP2
        else:
            JX = KX
            JY = KY
            for J in range(N):
                TEMP1 = ALPHA * X[JX]
                TEMP2 = 0
                Y[JY] += TEMP1 * A[J, J].real
                IX = JX
                IY = JY
                for I in range(J, N):
                    IX += INCX
                    IY += INCY
                    Y[IY] += TEMP1 * A[I, J]
                    TEMP2 += A[1, J].conjugate() * X[IX]
                Y[JY] += ALPHA * TEMP2
                JX += INCX
                JY += INCY
