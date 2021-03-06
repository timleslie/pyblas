# > \brief \b DGER
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def DGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
#
#       .. Scalar Arguments ..
#       DOUBLE PRECISION ALPHA
#       INTEGER INCX,INCY,LDA,M,N
#       ..
#       .. Array Arguments ..
#       DOUBLE PRECISION A(LDA,*),X(*),Y(*)
#       ..
#
#
# > \par Purpose:
#  =============
# >
# > \verbatim
# >
# > DGER   performs the rank 1 operation
# >
# >    A := alpha*x*y**T + A,
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
# >          ALPHA is DOUBLE PRECISION.
# >           On entry, ALPHA specifies the scalar alpha.
# > \endverbatim
# >
# > \param[in] X
# > \verbatim
# >          X is DOUBLE PRECISION array, dimension at least
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
# >          Y is DOUBLE PRECISION array, dimension at least
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
# >          A is DOUBLE PRECISION array, dimension ( LDA, N )
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
from xerbla import xerbla
from ..util import range_, slice_


def DGER(M, N, ALPHA, X, INCX, Y, INCY, A, LDA):
    #
    #  -- Reference BLAS level2 routine (version 3.7.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     December 2016
    #
    #     .. Scalar Arguments ..
    #   DOUBLE PRECISION ALPHA
    #   INTEGER INCX,INCY,LDA,M,N
    #     ..
    #     .. Array Arguments ..
    #   DOUBLE PRECISION A(LDA,*),X(*),Y(*)
    #     ..
    #
    #  =====================================================================

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
        xerbla("DGER  ", INFO)

    for J, JY in enumerate(range_(N, INCY)):
        A[:M, J] += ALPHA * X[slice_(M, INCX)] * Y[JY]
