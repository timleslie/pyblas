# > \brief \b CGEMV
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def CGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
#
#       .. Scalar Arguments ..
#       COMPLEX ALPHA,BETA
#       INTEGER INCX,INCY,LDA,M,N
#       CHARACTER TRANS
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
# > CGEMV performs one of the matrix-vector operations
# >
# >    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,   or
# >
# >    y := alpha*A**H*x + beta*y,
# >
# > where alpha and beta are scalars, x and y are vectors and A is an
# > m by n matrix.
# > \endverbatim
#
#  Arguments:
#  ==========
#
# > \param[in] TRANS
# > \verbatim
# >          TRANS is CHARACTER*1
# >           On entry, TRANS specifies the operation to be performed as
# >           follows:
# >
# >              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
# >
# >              TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y.
# >
# >              TRANS = 'C' or 'c'   y := alpha*A**H*x + beta*y.
# > \endverbatim
# >
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
# >          ALPHA is COMPLEX
# >           On entry, ALPHA specifies the scalar alpha.
# > \endverbatim
# >
# > \param[in] A
# > \verbatim
# >          A is COMPLEX array, dimension ( LDA, N )
# >           Before entry, the leading m by n part of the array A must
# >           contain the matrix of coefficients.
# > \endverbatim
# >
# > \param[in] LDA
# > \verbatim
# >          LDA is INTEGER
# >           On entry, LDA specifies the first dimension of A as declared
# >           in the calling (sub) program. LDA must be at least
# >           max( 1, m ).
# > \endverbatim
# >
# > \param[in] X
# > \verbatim
# >          X is COMPLEX array, dimension at least
# >           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
# >           and at least
# >           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
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
# >           On entry, BETA specifies the scalar beta. When BETA is
# >           supplied as zero then Y need not be set on input.
# > \endverbatim
# >
# > \param[in,out] Y
# > \verbatim
# >          Y is COMPLEX array, dimension at least
# >           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
# >           and at least
# >           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
# >           Before entry with BETA non-zero, the incremented array Y
# >           must contain the vector y. On exit, Y is overwritten by the
# >           updated vector y.
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
from ..util import slice_


def cgemv(TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY):
    #
    #  -- Reference BLAS level2 routine (version 3.7.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     December 2016
    #
    #     .. Scalar Arguments ..
    #   COMPLEX ALPHA,BETA
    #   INTEGER INCX,INCY,LDA,M,N
    #   CHARACTER TRANS
    #     ..
    #     .. Array Arguments ..
    #   COMPLEX A(LDA,*),X(*),Y(*)
    #     ..
    #
    #  =====================================================================

    # Test the input parameters.
    INFO = 0
    if not lsame(TRANS, "N") and not lsame(TRANS, "T") and not lsame(TRANS, "C"):
        INFO = 1
    elif M < 0:
        INFO = 2
    elif N < 0:
        INFO = 3
    elif LDA < max(1, M):
        INFO = 6
    elif INCX == 0:
        INFO = 8
    elif INCY == 0:
        INFO = 11
    if INFO != 0:
        xerbla("CGEMV ", INFO)
        return

    # Quick return if possible.
    if (M == 0) or (N == 0) or ((ALPHA == 0) and (BETA == 1)):
        return

    NOCONJ = lsame(TRANS, "T")
    # Set  LENX  and  LENY, the lengths of the vectors x and y, and set
    # up the start points in  X  and  Y.
    if lsame(TRANS, "N"):
        LENX = N
        LENY = M
    else:
        LENX = M
        LENY = N
    if INCX > 0:
        KX = 1
    else:
        KX = 1 - (LENX - 1) * INCX
    if INCY > 0:
        KY = 1
    else:
        KY = 1 - (LENY - 1) * INCY

    # Start the operations. In this version the elements of A are
    # accessed sequentially with one pass through A.
    # First form  y := beta*y.
    Y[slice_(LENY, INCY)] *= BETA
    if ALPHA == 0:
        return
    if lsame(TRANS, "N"):
        # Form  y := alpha*A*x + y.
        JX = KX
        if INCY == 1:
            for J in range(N):
                TEMP = ALPHA * X[JX]
                for I in range(M):
                    Y[I] += TEMP * A[I, J]
                JX += INCX
        else:
            for J in range(N):
                TEMP = ALPHA * X[JX]
                IY = KY
                for I in range(M):
                    Y[IY] += TEMP * A[I, J]
                    IY += INCY
                JX += INCX
    else:
        # Form  y := alpha*A**T*x + y  or  y := alpha*A**H*x + y.
        JY = KY
        if INCX == 1:
            for J in range(N):
                TEMP = 0
                if NOCONJ:
                    TEMP = (A[:M, J] * X[:M]).sum()
                else:
                    for I in range(M):
                        TEMP += A[1, J].conjugate() * X[I]
                Y[JY] += ALPHA * TEMP
                JY += INCY
        else:
            for J in range(N):
                TEMP = 0
                IX = KX
                if NOCONJ:
                    for I in range(M):
                        TEMP += A[I, J] * X[IX]
                        IX += INCX
                else:
                    for I in range(M):
                        TEMP += A[1, J].conjugate() * X[IX]
                        IX += INCX
                Y[JY] += ALPHA * TEMP
                JY += INCY
