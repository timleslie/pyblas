# > \brief \b DGBMV
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def DGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
#
#       .. Scalar Arguments ..
#       DOUBLE PRECISION ALPHA,BETA
#       INTEGER INCX,INCY,KL,KU,LDA,M,N
#       CHARACTER TRANS
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
# > DGBMV  performs one of the matrix-vector operations
# >
# >    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
# >
# > where alpha and beta are scalars, x and y are vectors and A is an
# > m by n band matrix, with kl sub-diagonals and ku super-diagonals.
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
# >              TRANS = 'C' or 'c'   y := alpha*A**T*x + beta*y.
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
# > \param[in] KL
# > \verbatim
# >          KL is INTEGER
# >           On entry, KL specifies the number of sub-diagonals of the
# >           matrix A. KL must satisfy  0 <= KL.
# > \endverbatim
# >
# > \param[in] KU
# > \verbatim
# >          KU is INTEGER
# >           On entry, KU specifies the number of super-diagonals of the
# >           matrix A. KU must satisfy  0 <= KU.
# > \endverbatim
# >
# > \param[in] ALPHA
# > \verbatim
# >          ALPHA is DOUBLE PRECISION.
# >           On entry, ALPHA specifies the scalar alpha.
# > \endverbatim
# >
# > \param[in] A
# > \verbatim
# >          A is DOUBLE PRECISION array, dimension ( LDA, N )
# >           Before entry, the leading ( kl + ku + 1 ) by n part of the
# >           array A must contain the matrix of coefficients, supplied
# >           column by column, with the leading diagonal of the matrix in
# >           row ( ku + 1 ) of the array, the first super-diagonal
# >           starting at position 2 in row ku, the first sub-diagonal
# >           starting at position 1 in row ( ku + 2 ), and so on.
# >           Elements in the array A that do not correspond to elements
# >           in the band matrix (such as the top left ku by ku triangle)
# >           are not referenced.
# >           The following program segment will transfer a band matrix
# >           from conventional full matrix storage to band storage:
# >
# >                 DO 20, J = 1, N
# >                    K = KU + 1 - J
# >                    DO 10, I = max( 1, J - KU ), min( M, J + KL )
# >                       A( K + I, J ) = matrix( I, J )
# >              10    CONTINUE
# >              20 CONTINUE
# > \endverbatim
# >
# > \param[in] LDA
# > \verbatim
# >          LDA is INTEGER
# >           On entry, LDA specifies the first dimension of A as declared
# >           in the calling (sub) program. LDA must be at least
# >           ( kl + ku + 1 ).
# > \endverbatim
# >
# > \param[in] X
# > \verbatim
# >          X is DOUBLE PRECISION array, dimension at least
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
# >          BETA is DOUBLE PRECISION.
# >           On entry, BETA specifies the scalar beta. When BETA is
# >           supplied as zero then Y need not be set on input.
# > \endverbatim
# >
# > \param[in,out] Y
# > \verbatim
# >          Y is DOUBLE PRECISION array, dimension at least
# >           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
# >           and at least
# >           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
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
# > \ingroup double_blas_level2
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


def DGBMV(TRANS, M, N, KL, KU, ALPHA, A, LDA, X, INCX, BETA, Y, INCY):
    #
    #  -- Reference BLAS level2 routine (version 3.7.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     December 2016
    #
    #     .. Scalar Arguments ..
    #   DOUBLE PRECISION ALPHA,BETA
    #   INTEGER INCX,INCY,KL,KU,LDA,M,N
    #   CHARACTER TRANS
    #     ..
    #     .. Array Arguments ..
    #   DOUBLE PRECISION A(LDA,*),X(*),Y(*)
    #     ..
    #
    #  =====================================================================
    #
    #     .. Parameters ..
    #   DOUBLE PRECISION ONE,ZERO
    #   PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
    #     ..
    #     .. Local Scalars ..
    #   DOUBLE PRECISION TEMP
    #   INTEGER I,INFO,IX,IY,J,JX,JY,K,KUP1,KX,KY,LENX,LENY
    #     ..
    #     .. External Functions ..
    #   LOGICAL LSAME
    #   EXTERNAL LSAME
    #     ..
    #     .. External Subroutines ..
    #   EXTERNAL XERBLA
    #     ..
    # .. Intrinsic Functions ..
    #   INTRINSIC MAX,MIN
    #     ..

    # Test the input parameters.
    INFO = 0
    if not lsame(TRANS, "N") and not lsame(TRANS, "T") and not lsame(TRANS, "C"):
        INFO = 1
    elif M < 0:
        INFO = 2
    elif N < 0:
        INFO = 3
    elif KL < 0:
        INFO = 4
    elif KU < 0:
        INFO = 5
    elif LDA < (KL + KU + 1):
        INFO = 8
    elif INCX == 0:
        INFO = 10
    elif INCY == 0:
        INFO = 13
    if INFO != 0:
        xerbla("DGBMV ", INFO)

    # Quick return if possible.
    if (M == 0) or (N == 0) or ((ALPHA == 0) and (BETA == 1)):
        return
    #
    #     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
    #     up the start points in  X  and  Y.
    #
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
    #
    #     Start the operations. In this version the elements of A are
    #     accessed sequentially with one pass through the band part of A.
    #
    # First form  y := beta*y.
    if INCY > 0:
        Y[: LENY * INCY : INCY] *= BETA
    else:
        Y[-(LENY - 1) * INCY :: INCY] *= BETA

    if ALPHA == 0:
        return
    KUP1 = KU + 1
    if lsame(TRANS, "N"):
        #
        #        Form  y := alpha*A*x + y.
        #
        JX = KX
        if INCY == 1:
            for J in range(N):
                TEMP = ALPHA * X[JX]
                K = KUP1 - J
                for I in range(max(1, J - KU) - 1, min(M, J + KL)):
                    Y[I] += TEMP * A[K + I, J]
                JX += INCX
        else:
            for J in range(N):
                TEMP = ALPHA * X[JX]
                IY = KY
                K = KUP1 - J
                for I in range(max(1, J - KU) - 1, min(M, J + KL)):
                    Y[IY] += TEMP * A[K + I, J]
                    IY += INCY
                JX += INCX
                if J > KU:
                    KY += INCY
    else:
        # Form  y := alpha*A**T*x + y.
        JY = KY
        if INCX == 1:
            for J in range(N):
                TEMP = 0
                K = KUP1 - J
                for I in range(max(1, J - KU) - 1, min(M, J + KL)):
                    TEMP += A[K + I, J] * X[I]
                Y[JY] += ALPHA * TEMP
                JY += INCY
        else:
            for J in range(N):
                TEMP = 0
                IX = KX
                K = KUP1 - J
                for I in range(max(1, J - KU) - 1, min(M, J + KL)):
                    TEMP += A[K + I, J] * X[IX]
                    IX += INCX
                Y[JY] += ALPHA * TEMP
                JY += INCY
                if J > KU:
                    KX += INCX
