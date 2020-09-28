# > \brief \b CTRSV
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def CTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
#
#       .. Scalar Arguments ..
#       INTEGER INCX,LDA,N
#       CHARACTER DIAG,TRANS,UPLO
#       ..
#       .. Array Arguments ..
#       COMPLEX A(LDA,*),X(*)
#       ..
#
#
# > \par Purpose:
#  =============
# >
# > \verbatim
# >
# > CTRSV  solves one of the systems of equations
# >
# >    A*x = b,   or   A**T*x = b,   or   A**H*x = b,
# >
# > where b and x are n element vectors and A is an n by n unit, or
# > non-unit, upper or lower triangular matrix.
# >
# > No test for singularity or near-singularity is included in this
# > routine. Such tests must be performed before calling this routine.
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
# >           On entry, TRANS specifies the equations to be solved as
# >           follows:
# >
# >              TRANS = 'N' or 'n'   A*x = b.
# >
# >              TRANS = 'T' or 't'   A**T*x = b.
# >
# >              TRANS = 'C' or 'c'   A**H*x = b.
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
# > \param[in] A
# > \verbatim
# >          A is COMPLEX array, dimension ( LDA, N )
# >           Before entry with  UPLO = 'U' or 'u', the leading n by n
# >           upper triangular part of the array A must contain the upper
# >           triangular matrix and the strictly lower triangular part of
# >           A is not referenced.
# >           Before entry with UPLO = 'L' or 'l', the leading n by n
# >           lower triangular part of the array A must contain the lower
# >           triangular matrix and the strictly upper triangular part of
# >           A is not referenced.
# >           Note that when  DIAG = 'U' or 'u', the diagonal elements of
# >           A are not referenced either, but are assumed to be unity.
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
# > \param[in,out] X
# > \verbatim
# >          X is COMPLEX array, dimension at least
# >           ( 1 + ( n - 1 )*abs( INCX ) ).
# >           Before entry, the incremented array X must contain the n
# >           element right-hand side vector b. On exit, X is overwritten
# >           with the solution vector x.
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
# > \ingroup complex_blas_level2
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


def CTRSV(UPLO, TRANS, DIAG, N, A, LDA, X, INCX):
    #
    #  -- Reference BLAS level2 routine (version 3.7.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     December 2016
    #
    #     .. Scalar Arguments ..
    #   INTEGER INCX,LDA,N
    #   CHARACTER DIAG,TRANS,UPLO
    #     ..
    #     .. Array Arguments ..
    #   COMPLEX A(LDA,*),X(*)
    #     ..
    #
    #  =====================================================================

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
    elif LDA < max(1, N):
        INFO = 6
    elif INCX == 0:
        INFO = 8
    if INFO != 0:
        xerbla("CTRSV ", INFO)
        return

    # Quick return if possible.
    if N == 0:
        return
    #
    NOCONJ = lsame(TRANS, "T")
    NOUNIT = lsame(DIAG, "N")
    #
    #     Set up the start point in X if the increment is not unity. This
    #     will be  ( N - 1 )*INCX  too small for descending loops.
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
        # Form  x := inv( A )*x.
        if lsame(UPLO, "U"):
            if INCX == 1:
                for J in range(N - 1, -1, -1):
                    if X[J] != 0:
                        if NOUNIT:
                            X[J] = X[J] / A[J, J]
                        TEMP = X[J]
                        for I in range(J - 2, -1, -1):
                            X[I] -= TEMP * A[I, J]
            else:
                JX = KX + (N - 1) * INCX
                for J in range(N - 1, -1, -1):
                    if X[JX] != 0:
                        if NOUNIT:
                            X[JX] = X[JX] / A[J, J]
                        TEMP = X[JX]
                        IX = JX
                        for I in range(J - 2, -1, -1):
                            IX -= INCX
                            X[IX] -= TEMP * A[I, J]
                    JX -= INCX
        else:
            if INCX == 1:
                for J in range(N):
                    if X[J] != 0:
                        if NOUNIT:
                            X[J] = X[J] / A[J, J]
                        TEMP = X[J]
                        for I in range(J, N):
                            X[I] -= TEMP * A[I, J]
            else:
                JX = KX
                for J in range(N):
                    if X[JX] != 0:
                        if NOUNIT:
                            X[JX] = X[JX] / A[J, J]
                        TEMP = X[JX]
                        IX = JX
                        for I in range(J, N):
                            IX += INCX
                            X[IX] -= TEMP * A[I, J]
                    JX += INCX
    else:
        # Form  x := inv( A**T )*x  or  x := inv( A**H )*x.
        if lsame(UPLO, "U"):
            if INCX == 1:
                for J in range(N):
                    TEMP = X[J]
                    if NOCONJ:
                        for I in range(J - 1):
                            TEMP -= A[I, J] * X[I]
                        if NOUNIT:
                            TEMP = TEMP / A[J, J]
                    else:
                        for I in range(J - 1):
                            TEMP -= A[1, J].conjugate() * X[I]
                        if NOUNIT:
                            TEMP = TEMP / (A[J, J]).conjugate()
                    X[J] = TEMP
            else:
                JX = KX
                for J in range(N):
                    IX = KX
                    TEMP = X[JX]
                    if NOCONJ:
                        for I in range(J - 1):
                            TEMP -= A[I, J] * X[IX]
                            IX += INCX
                        if NOUNIT:
                            TEMP = TEMP / A[J, J]
                    else:
                        for I in range(J - 1):
                            TEMP -= A[1, J].conjugate() * X[IX]
                            IX += INCX
                        if NOUNIT:
                            TEMP = TEMP / (A[J, J]).conjugate()
                    X[JX] = TEMP
                    JX += INCX
        else:
            if INCX == 1:
                for J in range(N - 1, -1, -1):
                    TEMP = X[J]
                    if NOCONJ:
                        for I in range(N - 1, J - 1, -1):
                            TEMP -= A[I, J] * X[I]
                        if NOUNIT:
                            TEMP = TEMP / A[J, J]
                    else:
                        for I in range(N - 1, J - 1, -1):
                            TEMP -= A[1, J].conjugate() * X[I]
                        if NOUNIT:
                            TEMP = TEMP / (A[J, J]).conjugate()
                    X[J] = TEMP
            else:
                KX += (N - 1) * INCX
                JX = KX
                for J in range(N - 1, -1, -1):
                    IX = KX
                    TEMP = X[JX]
                    if NOCONJ:
                        for I in range(N - 1, J - 1, -1):
                            TEMP -= A[I, J] * X[IX]
                            IX -= INCX
                        if NOUNIT:
                            TEMP = TEMP / A[J, J]
                    else:
                        for I in range(N - 1, J - 1, -1):
                            TEMP -= A[1, J].conjugate() * X[IX]
                            IX -= INCX
                        if NOUNIT:
                            TEMP = TEMP / (A[J, J]).conjugate()
                    X[JX] = TEMP
                    JX -= INCX
