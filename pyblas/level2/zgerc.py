# > \brief \b ZGERC
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def ZGERC(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
#
#       .. Scalar Arguments ..
#       COMPLEX*16 ALPHA
#       INTEGER INCX,INCY,LDA,M,N
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
# > ZGERC  performs the rank 1 operation
# >
# >    A := alpha*x*y**H + A,
# >
# > where alpha is a scalar, x is an m element vector, y is an n element
# > vector and A is an m by n matrix.
# > \endverbatim
#
#  Arguments:
#  ==========
#
# > \param[in] M
# > \verbatim
# >          M is INTEGER
# >           On entry, M specifies the number of rows of the matrix A.
# >           M must be at least zero.
# > \endverbatim
# >
# > \param[in] N
# > \verbatim
# >          N is INTEGER
# >           On entry, N specifies the number of columns of the matrix A.
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
# >           ( 1 + ( m - 1 )*abs( INCX ) ).
# >           Before entry, the incremented array X must contain the m
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
# >           Before entry, the leading m by n part of the array A must
# >           contain the matrix of coefficients. On exit, A is
# >           overwritten by the updated matrix.
# > \endverbatim
# >
# > \param[in] LDA
# > \verbatim
# >          LDA is INTEGER
# >           On entry, LDA specifies the first dimension of A as declared
# >           in the calling (sub) program. LDA must be at least
# >           max( 1, m ).
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
from xerbla import xerbla


def ZGERC(M, N, ALPHA, X, INCX, Y, INCY, A, LDA):
    #
    #  -- Reference BLAS level2 routine (version 3.7.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     December 2016
    #
    #     .. Scalar Arguments ..
    #   COMPLEX*16 ALPHA
    #   INTEGER INCX,INCY,LDA,M,N
    #     ..
    #     .. Array Arguments ..
    #   COMPLEX*16 A(LDA,*),X(*),Y(*)
    #     ..
    #
    #  =====================================================================
    #
    #     .. Parameters ..
    #   COMPLEX*16 ZERO
    #   PARAMETER (ZERO= (0.0D+0,0.0D+0))
    #     ..
    #     .. Local Scalars ..
    #   COMPLEX*16 TEMP
    #   INTEGER I,INFO,IX,J,JY,KX
    #     ..
    #     .. External Subroutines ..
    #   EXTERNAL XERBLA
    #     ..
    #     .. Intrinsic Functions ..
    #   INTRINSIC DCONJG,MAX
    #     ..

    # Test the input parameters.
    INFO = 0
    if M < 0:
        INFO = 1
    elif N < 0:
        INFO = 2
    elif INCX == 0:
        INFO = 5
    elif INCY == 0:
        INFO = 7
    elif LDA < max(1, M):
        INFO = 9
    if INFO != 0:
        xerbla("ZGERC ", INFO)
        return

    # Quick return if possible.
    if (M == 0) or (N == 0) or (ALPHA == 0):
        return
    #
    #     Start the operations. In this version the elements of A are
    #     accessed sequentially with one pass through A.
    #
    if INCY > 0:
        JY = 1
    else:
        JY = 1 - (N - 1) * INCY
    if INCX == 1:
        for J in range(N):
            if Y[JY] != 0:
                TEMP = ALPHA * Y[JY].conjugate()
                for I in range(M):
                    A[I, J] += X[I] * TEMP
            JY += INCY
    else:
        if INCX > 0:
            KX = 1
        else:
            KX = 1 - (M - 1) * INCX
        for J in range(N):
            if Y[JY] != 0:
                TEMP = ALPHA * Y[JY].conjugate()
                IX = KX
                for I in range(M):
                    A[I, J] += X[IX] * TEMP
                    IX += INCX
            JY += INCY
