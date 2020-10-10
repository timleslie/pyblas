# > \brief \b ZHER2
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def ZHER2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)
#
#       .. Scalar Arguments ..
#       COMPLEX*16 ALPHA
#       INTEGER INCX,INCY,LDA,N
#       CHARACTER UPLO
#       ..
#       .. Array Arguments ..
#       COMPLEX*16 A(LDA,*),X(*),Y(*)
#       ..
#
#
# > \par Purpose:
#  =============
# >
# > \verbatim
# >
# > ZHER2  performs the hermitian rank 2 operation
# >
# >    A := alpha*x*y**H + conjg( alpha )*y*x**H + A,
# >
# > where alpha is a scalar, x and y are n element vectors and A is an n
# > by n hermitian matrix.
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
# >          ALPHA is COMPLEX*16
# >           On entry, ALPHA specifies the scalar alpha.
# > \endverbatim
# >
# > \param[in] X
# > \verbatim
# >          X is COMPLEX*16 array, dimension at least
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
# >          Y is COMPLEX*16 array, dimension at least
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
# >          A is COMPLEX*16 array, dimension ( LDA, N )
# >           Before entry with  UPLO = 'U' or 'u', the leading n by n
# >           upper triangular part of the array A must contain the upper
# >           triangular part of the hermitian matrix and the strictly
# >           lower triangular part of A is not referenced. On exit, the
# >           upper triangular part of the array A is overwritten by the
# >           upper triangular part of the updated matrix.
# >           Before entry with UPLO = 'L' or 'l', the leading n by n
# >           lower triangular part of the array A must contain the lower
# >           triangular part of the hermitian matrix and the strictly
# >           upper triangular part of A is not referenced. On exit, the
# >           lower triangular part of the array A is overwritten by the
# >           lower triangular part of the updated matrix.
# >           Note that the imaginary parts of the diagonal elements need
# >           not be set, they are assumed to be zero, and on exit they
# >           are set to zero.
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
# > \ingroup complex16_blas_level2
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


def ZHER2(UPLO, N, ALPHA, X, INCX, Y, INCY, A, LDA):
    #
    #  -- Reference BLAS level2 routine (version 3.7.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     December 2016
    #
    #     .. Scalar Arguments ..
    #   COMPLEX*16 ALPHA
    #   INTEGER INCX,INCY,LDA,N
    #   CHARACTER UPLO
    #     ..
    #     .. Array Arguments ..
    #   COMPLEX*16 A(LDA,*),X(*),Y(*)
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
        xerbla("ZHER2 ", INFO)
        return

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
                if (X[J] != 0) or (Y[J] != 0):
                    TEMP1 = ALPHA * (Y[J]).conjugate()
                    TEMP2 = (ALPHA * X[J]).conjugate()
                    for I in range(J - 1):
                        A[I, J] += X[I] * TEMP1 + Y[I] * TEMP2
                    A[J, J] = A[J, J].real + (X[J] * TEMP1 + Y[J] * TEMP2).real
                else:
                    A[J, J] = A[J, J].real
        else:
            for J in range(N):
                if (X[JX] != 0) or (Y[JY] != 0):
                    TEMP1 = ALPHA * Y[JY].conjugate()
                    TEMP2 = (ALPHA * X[JX]).conjugate()
                    IX = KX
                    IY = KY
                    for I in range(J - 1):
                        A[I, J] += X[IX] * TEMP1 + Y[IY] * TEMP2
                        IX += INCX
                        IY += INCY
                    A[J, J] = A[J, J].real + (X[JX] * TEMP1 + Y[JY] * TEMP2).real
                else:
                    A[J, J] = A[J, J].real
                JX += INCX
                JY += INCY
    else:
        # Form  A  when A is stored in the lower triangle.
        if (INCX == 1) and (INCY == 1):
            for J in range(N):
                if (X[J] != 0) or (Y[J] != 0):
                    TEMP1 = ALPHA * (Y[J]).conjugate()
                    TEMP2 = (ALPHA * X[J]).conjugate()
                    A[J, J] = A[J, J].real + (X[J] * TEMP1 + Y[J] * TEMP2).real
                    for I in range(J, N):
                        A[I, J] += X[I] * TEMP1 + Y[I] * TEMP2
                else:
                    A[J, J] = A[J, J].real
        else:
            for J in range(N):
                if (X[JX] != 0) or (Y[JY] != 0):
                    TEMP1 = ALPHA * Y[JY].conjugate()
                    TEMP2 = (ALPHA * X[JX]).conjugate()
                    A[J, J] = A[J, J].real + (X[JX] * TEMP1 + Y[JY] * TEMP2).real
                    IX = JX
                    IY = JY
                    for I in range(J, N):
                        IX += INCX
                        IY += INCY
                        A[I, J] += X[IX] * TEMP1 + Y[IY] * TEMP2
                else:
                    A[J, J] = A[J, J].real
                JX += INCX
                JY += INCY
