# > \brief \b CSYRK
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def CSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
#
#       .. Scalar Arguments ..
#       COMPLEX ALPHA,BETA
#       INTEGER K,LDA,LDC,N
#       CHARACTER TRANS,UPLO
#       ..
#       .. Array Arguments ..
#       COMPLEX A(LDA,*),C(LDC,*)
#       ..
#
#
# > \par Purpose:
#  =============
# >
# > \verbatim
# >
# > CSYRK  performs one of the symmetric rank k operations
# >
# >    C := alpha*A*A**T + beta*C,
# >
# > or
# >
# >    C := alpha*A**T*A + beta*C,
# >
# > where  alpha and beta  are scalars,  C is an  n by n symmetric matrix
# > and  A  is an  n by k  matrix in the first case and a  k by n  matrix
# > in the second case.
# > \endverbatim
#
#  Arguments:
#  ==========
#
# > \param[in] UPLO
# > \verbatim
# >          UPLO is CHARACTER*1
# >           On  entry,   UPLO  specifies  whether  the  upper  or  lower
# >           triangular  part  of the  array  C  is to be  referenced  as
# >           follows:
# >
# >              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
# >                                  is to be referenced.
# >
# >              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
# >                                  is to be referenced.
# > \endverbatim
# >
# > \param[in] TRANS
# > \verbatim
# >          TRANS is CHARACTER*1
# >           On entry,  TRANS  specifies the operation to be performed as
# >           follows:
# >
# >              TRANS = 'N' or 'n'   C := alpha*A*A**T + beta*C.
# >
# >              TRANS = 'T' or 't'   C := alpha*A**T*A + beta*C.
# > \endverbatim
# >
# > \param[in] N
# > \verbatim
# >          N is INTEGER
# >           On entry,  N specifies the order of the matrix C.  N must be
# >           at least zero.
# > \endverbatim
# >
# > \param[in] K
# > \verbatim
# >          K is INTEGER
# >           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
# >           of  columns   of  the   matrix   A,   and  on   entry   with
# >           TRANS = 'T' or 't',  K  specifies  the number of rows of the
# >           matrix A.  K must be at least zero.
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
# >          A is COMPLEX array, dimension ( LDA, ka ), where ka is
# >           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
# >           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
# >           part of the array  A  must contain the matrix  A,  otherwise
# >           the leading  k by n  part of the array  A  must contain  the
# >           matrix A.
# > \endverbatim
# >
# > \param[in] LDA
# > \verbatim
# >          LDA is INTEGER
# >           On entry, LDA specifies the first dimension of A as declared
# >           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
# >           then  LDA must be at least  max( 1, n ), otherwise  LDA must
# >           be at least  max( 1, k ).
# > \endverbatim
# >
# > \param[in] BETA
# > \verbatim
# >          BETA is COMPLEX
# >           On entry, BETA specifies the scalar beta.
# > \endverbatim
# >
# > \param[in,out] C
# > \verbatim
# >          C is COMPLEX array, dimension ( LDC, N )
# >           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
# >           upper triangular part of the array C must contain the upper
# >           triangular part  of the  symmetric matrix  and the strictly
# >           lower triangular part of C is not referenced.  On exit, the
# >           upper triangular part of the array  C is overwritten by the
# >           upper triangular part of the updated matrix.
# >           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
# >           lower triangular part of the array C must contain the lower
# >           triangular part  of the  symmetric matrix  and the strictly
# >           upper triangular part of C is not referenced.  On exit, the
# >           lower triangular part of the array  C is overwritten by the
# >           lower triangular part of the updated matrix.
# > \endverbatim
# >
# > \param[in] LDC
# > \verbatim
# >          LDC is INTEGER
# >           On entry, LDC specifies the first dimension of C as declared
# >           in  the  calling  (sub)  program.   LDC  must  be  at  least
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
# > \ingroup complex_blas_level3
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


def CSYRK(UPLO, TRANS, N, K, ALPHA, A, LDA, BETA, C, LDC):
    #
    #  -- Reference BLAS level3 routine (version 3.7.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     December 2016
    #
    #     .. Scalar Arguments ..
    #   COMPLEX ALPHA,BETA
    #   INTEGER K,LDA,LDC,N
    #   CHARACTER TRANS,UPLO
    #     ..
    #     .. Array Arguments ..
    #   COMPLEX A(LDA,*),C(LDC,*)
    #     ..
    #
    #  =====================================================================
    #
    #     .. External Functions ..
    #   LOGICAL LSAME
    #   EXTERNAL LSAME
    #     ..
    #     .. External Subroutines ..
    #   EXTERNAL XERBLA
    #     ..
    #     .. Intrinsic Functions ..
    #   INTRINSIC MAX
    #     ..
    #     .. Local Scalars ..
    #   COMPLEX TEMP
    #   INTEGER I,INFO,J,L,NROWA
    #   LOGICAL UPPER
    #     ..
    #     .. Parameters ..
    #   COMPLEX ONE
    #   PARAMETER (ONE= (1.0E+0,0.0E+0))
    #   COMPLEX ZERO
    #   PARAMETER (ZERO= (0.0E+0,0.0E+0))
    #     ..

    # Test the input parameters.
    if lsame(TRANS, "N"):
        NROWA = N
    else:
        NROWA = K
    UPPER = lsame(UPLO, "U")
    #
    INFO = 0
    if (not UPPER) and (not lsame(UPLO, "L")):
        INFO = 1
    elif (not lsame(TRANS, "N")) and (not lsame(TRANS, "T")):
        INFO = 2
    elif N < 0:
        INFO = 3
    elif K < 0:
        INFO = 4
    elif LDA < max(1, NROWA):
        INFO = 7
    elif LDC < max(1, N):
        INFO = 10
    if INFO != 0:
        xerbla("CSYRK ", INFO)
        return

    # Quick return if possible.
    if (N == 0) or (((ALPHA == 0) or (K == 0)) and (BETA == 1)):
        return

    # And when  alpha==zero.
    if ALPHA == 0:
        if UPPER:
            if BETA == 0:
                for J in range(N):
                    for I in range(J):
                        C[I, J] = 0
            else:
                for J in range(N):
                    for I in range(J):
                        C[I, J] *= BETA
        else:
            if BETA == 0:
                for J in range(N):
                    for I in range(J - 1, N):
                        C[I, J] = 0
            else:
                for J in range(N):
                    for I in range(J - 1, N):
                        C[I, J] *= BETA
        return

    # Start the operations.
    if lsame(TRANS, "N"):
        #
        #        Form  C := alpha*A*A**T + beta*C.
        #
        if UPPER:
            for J in range(N):
                if BETA == 0:
                    for I in range(J):
                        C[I, J] = 0
                elif BETA != 1:
                    for I in range(J):
                        C[I, J] *= BETA
                for L in range(K):
                    if A[J, L] != 0:
                        TEMP = ALPHA * A[J, L]
                        for I in range(J):
                            C[I, J] += TEMP * A[I, L]
        else:
            for J in range(N):
                if BETA == 0:
                    for I in range(J - 1, N):
                        C[I, J] = 0
                elif BETA != 1:
                    for I in range(J - 1, N):
                        C[I, J] *= BETA
                for L in range(K):
                    if A[J, L] != 0:
                        TEMP = ALPHA * A[J, L]
                        for I in range(J - 1, N):
                            C[I, J] += TEMP * A[I, L]
    else:
        # Form  C := alpha*A**T*A + beta*C.
        if UPPER:
            for J in range(N):
                for I in range(J):
                    TEMP = 0
                    for L in range(K):
                        TEMP += A[L, I] * A[L, J]
                    if BETA == 0:
                        C[I, J] = ALPHA * TEMP
                    else:
                        C[I, J] = ALPHA * TEMP + BETA * C[I, J]
        else:
            for J in range(N):
                for I in range(J - 1, N):
                    TEMP = 0
                    for L in range(K):
                        TEMP += A[L, I] * A[L, J]
                    if BETA == 0:
                        C[I, J] = ALPHA * TEMP
                    else:
                        C[I, J] = ALPHA * TEMP + BETA * C[I, J]
