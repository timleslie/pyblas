# > \brief \b SSYR2
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def SSYR2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)
#
#       .. Scalar Arguments ..
#       REAL ALPHA
#       INTEGER INCX,INCY,LDA,N
#       CHARACTER UPLO
#       ..
#       .. Array Arguments ..
#       REAL A(LDA,*),X(*),Y(*)
#       ..
#
#
# > \par Purpose:
#  =============
# >
# > \verbatim
# >
# > SSYR2  performs the symmetric rank 2 operation
# >
# >    A := alpha*x*y**T + alpha*y*x**T + A,
# >
# > where alpha is a scalar, x and y are n element vectors and A is an n
# > by n symmetric matrix.
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
# > \param[in] Y
# > \verbatim
# >          Y is REAL array, dimension at least
# >           ( 1 + ( n - 1 )*abs( INCY ) ).
# >           Before entry, the incremented array Y must contain the n
# >           element vector y.
# > \endverbatim
# >
# > \param[in] INCY
# > \verbatim
# >          INCY is INTEGER
# >           On entry, INCY specifies the increment for the elements of
# >           Y. INCY must not be zero.
# > \endverbatim
# >
# > \param[in,out] A
# > \verbatim
# >          A is REAL array, dimension ( LDA, N )
# >           Before entry with  UPLO = 'U' or 'u', the leading n by n
# >           upper triangular part of the array A must contain the upper
# >           triangular part of the symmetric matrix and the strictly
# >           lower triangular part of A is not referenced. On exit, the
# >           upper triangular part of the array A is overwritten by the
# >           upper triangular part of the updated matrix.
# >           Before entry with UPLO = 'L' or 'l', the leading n by n
# >           lower triangular part of the array A must contain the lower
# >           triangular part of the symmetric matrix and the strictly
# >           upper triangular part of A is not referenced. On exit, the
# >           lower triangular part of the array A is overwritten by the
# >           lower triangular part of the updated matrix.
# > \endverbatim
# >
# > \param[in] LDA
# > \verbatim
# >          LDA is INTEGER
# >           On entry, LDA specifies the first dimension of A as declared
# >           in the calling (sub) program. LDA must be at least
# >           max( 1, n ).
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


def SSYR2(UPLO, N, ALPHA, X, INCX, Y, INCY, A, LDA):
    #
    #  -- Reference BLAS level2 routine (version 3.7.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     December 2016
    #
    #     .. Scalar Arguments ..
    #   REAL ALPHA
    #   INTEGER INCX,INCY,LDA,N
    #   CHARACTER UPLO
    #     ..
    #     .. Array Arguments ..
    #   REAL A(LDA,*),X(*),Y(*)
    #     ..
    #
    #  =====================================================================

    # Test the input parameters.
    INFO = 0
    if not lsame(UPLO, "U") and not lsame(UPLO, "L"):
        INFO = 1
    elif N < 0:
        INFO = 2
    elif INCX == 0:
        INFO = 5
    elif INCY == 0:
        INFO = 7
    elif LDA < max(1, N):
        INFO = 9
    if INFO != 0:
        xerbla("SSYR2 ", INFO)

    # Quick return if possible.
    if (N == 0) or (ALPHA == 0):
        return
    #
    #     Set up the start points in X and Y if the increments are not both
    #     unity.
    #
    if (INCX != 1) or (INCY != 1):
        if INCX > 0:
            KX = 1
        else:
            KX = 1 - (N - 1) * INCX
        if INCY > 0:
            KY = 1
        else:
            KY = 1 - (N - 1) * INCY
        JX = KX
        JY = KY
    #
    #     Start the operations. In this version the elements of A are
    #     accessed sequentially with one pass through the triangular part
    #     of A.
    #
    if lsame(UPLO, "U"):
        #
        #        Form  A  when A is stored in the upper triangle.
        #
        if (INCX == 1) and (INCY == 1):
            for J in range(N):
                TEMP1 = ALPHA * Y[J]
                TEMP2 = ALPHA * X[J]
                for I in range(J):
                    A[I, J] += X[I] * TEMP1 + Y[I] * TEMP2
        else:
            for J in range(N):
                if (X[JX] != 0) or (Y[JY] != 0):
                    TEMP1 = ALPHA * Y[JY]
                    TEMP2 = ALPHA * X[JX]
                    IX = KX
                    IY = KY
                    for I in range(J):
                        A[I, J] += X[IX] * TEMP1 + Y[IY] * TEMP2
                        IX += INCX
                        IY += INCY
                JX += INCX
                JY += INCY
    else:
        # Form  A  when A is stored in the lower triangle.
        if (INCX == 1) and (INCY == 1):
            for J in range(N):
                if (X[J] != 0) or (Y[J] != 0):
                    TEMP1 = ALPHA * Y[J]
                    TEMP2 = ALPHA * X[J]
                    for I in range(J - 1, N):
                        A[I, J] += X[I] * TEMP1 + Y[I] * TEMP2
        else:
            for J in range(N):
                if (X[JX] != 0) or (Y[JY] != 0):
                    TEMP1 = ALPHA * Y[JY]
                    TEMP2 = ALPHA * X[JX]
                    IX = JX
                    IY = JY
                    for I in range(J - 1, N):
                        A[I, J] += X[IX] * TEMP1 + Y[IY] * TEMP2
                        IX += INCX
                        IY += INCY
                JX += INCX
                JY += INCY
