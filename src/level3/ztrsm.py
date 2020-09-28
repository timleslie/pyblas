# > \brief \b ZTRSM
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
#
#       .. Scalar Arguments ..
#       COMPLEX*16 ALPHA
#       INTEGER LDA,LDB,M,N
#       CHARACTER DIAG,SIDE,TRANSA,UPLO
#       ..
#       .. Array Arguments ..
#       COMPLEX*16 A(LDA,*),B(LDB,*)
#       ..
#
#
# > \par Purpose:
#  =============
# >
# > \verbatim
# >
# > ZTRSM  solves one of the matrix equations
# >
# >    op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
# >
# > where alpha is a scalar, X and B are m by n matrices, A is a unit, or
# > non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
# >
# >    op( A ) = A   or   op( A ) = A**T   or   op( A ) = A**H.
# >
# > The matrix X is overwritten on B.
# > \endverbatim
#
#  Arguments:
#  ==========
#
# > \param[in] SIDE
# > \verbatim
# >          SIDE is CHARACTER*1
# >           On entry, SIDE specifies whether op( A ) appears on the left
# >           or right of X as follows:
# >
# >              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
# >
# >              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
# > \endverbatim
# >
# > \param[in] UPLO
# > \verbatim
# >          UPLO is CHARACTER*1
# >           On entry, UPLO specifies whether the matrix A is an upper or
# >           lower triangular matrix as follows:
# >
# >              UPLO = 'U' or 'u'   A is an upper triangular matrix.
# >
# >              UPLO = 'L' or 'l'   A is a lower triangular matrix.
# > \endverbatim
# >
# > \param[in] TRANSA
# > \verbatim
# >          TRANSA is CHARACTER*1
# >           On entry, TRANSA specifies the form of op( A ) to be used in
# >           the matrix multiplication as follows:
# >
# >              TRANSA = 'N' or 'n'   op( A ) = A.
# >
# >              TRANSA = 'T' or 't'   op( A ) = A**T.
# >
# >              TRANSA = 'C' or 'c'   op( A ) = A**H.
# > \endverbatim
# >
# > \param[in] DIAG
# > \verbatim
# >          DIAG is CHARACTER*1
# >           On entry, DIAG specifies whether or not A is unit triangular
# >           as follows:
# >
# >              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
# >
# >              DIAG = 'N' or 'n'   A is not assumed to be unit
# >                                  triangular.
# > \endverbatim
# >
# > \param[in] M
# > \verbatim
# >          M is INTEGER
# >           On entry, M specifies the number of rows of B. M must be at
# >           least zero.
# > \endverbatim
# >
# > \param[in] N
# > \verbatim
# >          N is INTEGER
# >           On entry, N specifies the number of columns of B.  N must be
# >           at least zero.
# > \endverbatim
# >
# > \param[in] ALPHA
# > \verbatim
# >          ALPHA is COMPLEX*16
# >           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
# >           zero then  A is not referenced and  B need not be set before
# >           entry.
# > \endverbatim
# >
# > \param[in] A
# > \verbatim
# >          A is COMPLEX*16 array, dimension ( LDA, k ),
# >           where k is m when SIDE = 'L' or 'l'
# >             and k is n when SIDE = 'R' or 'r'.
# >           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
# >           upper triangular part of the array  A must contain the upper
# >           triangular matrix  and the strictly lower triangular part of
# >           A is not referenced.
# >           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
# >           lower triangular part of the array  A must contain the lower
# >           triangular matrix  and the strictly upper triangular part of
# >           A is not referenced.
# >           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
# >           A  are not referenced either,  but are assumed to be  unity.
# > \endverbatim
# >
# > \param[in] LDA
# > \verbatim
# >          LDA is INTEGER
# >           On entry, LDA specifies the first dimension of A as declared
# >           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
# >           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
# >           then LDA must be at least max( 1, n ).
# > \endverbatim
# >
# > \param[in,out] B
# > \verbatim
# >          B is COMPLEX*16 array, dimension ( LDB, N )
# >           Before entry,  the leading  m by n part of the array  B must
# >           contain  the  right-hand  side  matrix  B,  and  on exit  is
# >           overwritten by the solution matrix  X.
# > \endverbatim
# >
# > \param[in] LDB
# > \verbatim
# >          LDB is INTEGER
# >           On entry, LDB specifies the first dimension of B as declared
# >           in  the  calling  (sub)  program.   LDB  must  be  at  least
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
# > \ingroup complex16_blas_level3
#
# > \par Further Details:
#  =====================
# >
# > \verbatim
# >
# >  Level 3 Blas routine.
# >
# >  -- Written on 8-February-1989.
# >     Jack Dongarra, Argonne National Laboratory.
# >     Iain Duff, AERE Harwell.
# >     Jeremy Du Croz, Numerical Algorithms Group Ltd.
# >     Sven Hammarling, Numerical Algorithms Group Ltd.
# > \endverbatim
# >
#  =====================================================================
from util import lsame
from xerbla import xerbla


def ZTRSM(SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB):
    #
    #  -- Reference BLAS level3 routine (version 3.7.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     December 2016
    #
    #     .. Scalar Arguments ..
    #   COMPLEX*16 ALPHA
    #   INTEGER LDA,LDB,M,N
    #   CHARACTER DIAG,SIDE,TRANSA,UPLO
    #     ..
    #     .. Array Arguments ..
    #   COMPLEX*16 A(LDA,*),B(LDB,*)
    #     ..
    #
    #  =====================================================================

    # Test the input parameters.
    LSIDE = lsame(SIDE, "L")
    if LSIDE:
        NROWA = M
    else:
        NROWA = N
    NOCONJ = lsame(TRANSA, "T")
    NOUNIT = lsame(DIAG, "N")
    UPPER = lsame(UPLO, "U")

    INFO = 0
    if (not LSIDE) and (not lsame(SIDE, "R")):
        INFO = 1
    elif (not UPPER) and (not lsame(UPLO, "L")):
        INFO = 2
    elif (
        (not lsame(TRANSA, "N"))
        and (not lsame(TRANSA, "T"))
        and (not lsame(TRANSA, "C"))
    ):
        INFO = 3
    elif (not lsame(DIAG, "U")) and (not lsame(DIAG, "N")):
        INFO = 4
    elif M < 0:
        INFO = 5
    elif N < 0:
        INFO = 6
    elif LDA < max(1, NROWA):
        INFO = 9
    elif LDB < max(1, M):
        INFO = 11
    if INFO != 0:
        xerbla("ZTRSM ", INFO)
        return

    # Quick return if possible.
    if M == 0 or N == 0:
        return

    # And when  alpha==zero.
    if ALPHA == 0:
        for J in range(N):
            for I in range(M):
                B[I, J] = 0
        return

    # Start the operations.
    if LSIDE:
        if lsame(TRANSA, "N"):
            #
            #           Form  B := alpha*inv( A )*B.
            #
            if UPPER:
                for J in range(N):
                    if ALPHA != 1:
                        for I in range(M):
                            B[I, J] = ALPHA * B[I, J]
                    for K in range(M - 1, -1, -1):
                        if B[K, J] != 0:
                            if NOUNIT:
                                B[K, J] = B[K, J] / A[K, K]
                            for I in range(K - 1):
                                B[I, J] = B[I, J] - B[K, J] * A[I, K]
            else:
                for J in range(N):
                    if ALPHA != 1:
                        for I in range(M):
                            B[I, J] = ALPHA * B[I, J]
                    for K in range(M):
                        if B[K, J] != 0:
                            if NOUNIT:
                                B[K, J] = B[K, J] / A[K, K]
                            for I in range(K, M):
                                B[I, J] = B[I, J] - B[K, J] * A[I, K]
        else:
            # Form  B := alpha*inv( A**T )*B
            # or    B := alpha*inv( A**H )*B.
            if UPPER:
                for J in range(N):
                    for I in range(M):
                        TEMP = ALPHA * B[I, J]
                        if NOCONJ:
                            for K in range(I - 1):
                                TEMP -= A[K, I] * B[K, J]
                            if NOUNIT:
                                TEMP = TEMP / A[I, I]
                        else:
                            for K in range(I - 1):
                                TEMP -= A[K, I].conjugate() * B[K, J]
                            if NOUNIT:
                                TEMP = TEMP / A[I, I].conjugate()
                        B[I, J] = TEMP
            else:
                for J in range(N):
                    for I in range(M - 1, -1, -1):
                        TEMP = ALPHA * B[I, J]
                        if NOCONJ:
                            for K in range(I, M):
                                TEMP -= A[K, I] * B[K, J]
                            if NOUNIT:
                                TEMP = TEMP / A[I, I]
                        else:
                            for K in range(I, M):
                                TEMP -= A[K, I].conjugate() * B[K, J]
                            if NOUNIT:
                                TEMP = TEMP / A[I, I].conjugate()
                        B[I, J] = TEMP
    else:
        if lsame(TRANSA, "N"):
            # Form  B := alpha*B*inv( A ).
            if UPPER:
                for J in range(N):
                    if ALPHA != 1:
                        for I in range(M):
                            B[I, J] = ALPHA * B[I, J]
                    for K in range(J - 1):
                        if A[K, J] != 0:
                            for I in range(M):
                                B[I, J] = B[I, J] - A[K, J] * B[I, K]
                    if NOUNIT:
                        TEMP = 1 / A[J, J]
                        for I in range(M):
                            B[I, J] = TEMP * B[I, J]
            else:
                for J in range(N - 1, -1, -1):
                    if ALPHA != 1:
                        for I in range(M):
                            B[I, J] = ALPHA * B[I, J]
                    for K in range(J, N):
                        if A[K, J] != 0:
                            for I in range(M):
                                B[I, J] = B[I, J] - A[K, J] * B[I, K]
                    if NOUNIT:
                        TEMP = 1 / A[J, J]
                        for I in range(M):
                            B[I, J] = TEMP * B[I, J]
        else:
            # Form  B := alpha*B*inv( A**T )
            # or    B := alpha*B*inv( A**H ).
            if UPPER:
                for K in range(N - 1, -1, -1):
                    if NOUNIT:
                        if NOCONJ:
                            TEMP = 1 / A[K, K]
                        else:
                            TEMP = 1 / A[K, K].conjugate()
                        for I in range(M):
                            B[I, K] = TEMP * B[I, K]
                    for K in range(K - 1):
                        if A[J, K] != 0:
                            if NOCONJ:
                                TEMP = A[J, K]
                            else:
                                TEMP = A[J, K].conjugate()
                            for I in range(M):
                                B[I, J] = B[I, J] - TEMP * B[I, K]
                    if ALPHA != 1:
                        for I in range(M):
                            B[I, K] = ALPHA * B[I, K]
            else:
                for K in range(N):
                    if NOUNIT:
                        if NOCONJ:
                            TEMP = 1 / A[K, K]
                        else:
                            TEMP = 1 / A[K, K].conjugate()
                        for I in range(M):
                            B[I, K] = TEMP * B[I, K]
                    for J in range(K, N):
                        if A[J, K] != 0:
                            if NOCONJ:
                                TEMP = A[J, K]
                            else:
                                TEMP = A[J, K].conjugate()
                            for I in range(M):
                                B[I, J] = B[I, J] - TEMP * B[I, K]
                    if ALPHA != 1:
                        for I in range(M):
                            B[I, K] = ALPHA * B[I, K]
