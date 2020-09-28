# > \brief \b ZHER
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def ZHER(UPLO,N,ALPHA,X,INCX,A,LDA)
#
#       .. Scalar Arguments ..
#       DOUBLE PRECISION ALPHA
#       INTEGER INCX,LDA,N
#       CHARACTER UPLO
#       ..
#       .. Array Arguments ..
#       COMPLEX*16 A(LDA,*),X(*)
#       ..
#
#
# > \par Purpose:
#  =============
# >
# > \verbatim
# >
# > ZHER   performs the hermitian rank 1 operation
# >
# >    A := alpha*x*x**H + A,
# >
# > where alpha is a real scalar, x is an n element vector and A is an
# > n by n hermitian matrix.
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
# >          ALPHA is DOUBLE PRECISION.
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


def ZHER(UPLO, N, ALPHA, X, INCX, A, LDA):
    #
    #  -- Reference BLAS level2 routine (version 3.7.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     December 2016
    #
    #     .. Scalar Arguments ..
    #   DOUBLE PRECISION ALPHA
    #   INTEGER INCX,LDA,N
    #   CHARACTER UPLO
    #     ..
    #     .. Array Arguments ..
    #   COMPLEX*16 A(LDA,*),X(*)
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
    elif LDA < max(1, N):
        INFO = 7
    if INFO != 0:
        xerbla("ZHER  ", INFO)

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
    #     Start the operations. In this version the elements of A are
    #     accessed sequentially with one pass through the triangular part
    #     of A.
    #
    if lsame(UPLO, "U"):
        #
        #        Form  A  when A is stored in upper triangle.
        #
        if INCX == 1:
            for J in range(N):
                if X[J] != 0:
                    TEMP = ALPHA * (X[J]).conjugate()
                    for I in range(J - 1):
                        A[I, J] += X[I] * TEMP
                    A[J, J] = A[J, J].real + (X[J] * TEMP).real
                else:
                    A[J, J] = A[J, J].real
        else:
            JX = KX
            for J in range(N):
                if X[JX] != 0:
                    TEMP = ALPHA * (X[JX]).conjugate()
                    IX = KX
                    for I in range(J - 1):
                        A[I, J] += X[IX] * TEMP
                        IX += INCX
                    A[J, J] = A[J, J].real + (X[JX] * TEMP).real
                else:
                    A[J, J] = A[J, J].real
                JX += INCX
    else:
        # Form  A  when A is stored in lower triangle.
        if INCX == 1:
            for J in range(N):
                if X[J] != 0:
                    TEMP = ALPHA * (X[J]).conjugate()
                    A[J, J] = A[J, J].real + (TEMP * X[J]).real
                    for I in range(J, N):
                        A[I, J] += X[I] * TEMP
                else:
                    A[J, J] = A[J, J].real
        else:
            JX = KX
            for J in range(N):
                if X[JX] != 0:
                    TEMP = ALPHA * (X[JX]).conjugate()
                    A[J, J] = A[J, J].real + (TEMP * X[JX]).real
                    IX = JX
                    for I in range(J, N):
                        IX += INCX
                        A[I, J] += X[IX] * TEMP
                else:
                    A[J, J] = A[J, J].real
                JX += INCX
