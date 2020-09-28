# > \brief \b CHER2K
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def CHER2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
#
#       .. Scalar Arguments ..
#       COMPLEX ALPHA
#       REAL BETA
#       INTEGER K,LDA,LDB,LDC,N
#       CHARACTER TRANS,UPLO
#       ..
#       .. Array Arguments ..
#       COMPLEX A(LDA,*),B(LDB,*),C(LDC,*)
#       ..
#
#
# > \par Purpose:
#  =============
# >
# > \verbatim
# >
# > CHER2K  performs one of the hermitian rank 2k operations
# >
# >    C := alpha*A*B**H + conjg( alpha )*B*A**H + beta*C,
# >
# > or
# >
# >    C := alpha*A**H*B + conjg( alpha )*B**H*A + beta*C,
# >
# > where  alpha and beta  are scalars with  beta  real,  C is an  n by n
# > hermitian matrix and  A and B  are  n by k matrices in the first case
# > and  k by n  matrices in the second case.
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
# >              TRANS = 'N' or 'n'    C := alpha*A*B**H          +
# >                                         conjg( alpha )*B*A**H +
# >                                         beta*C.
# >
# >              TRANS = 'C' or 'c'    C := alpha*A**H*B          +
# >                                         conjg( alpha )*B**H*A +
# >                                         beta*C.
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
# >           of  columns  of the  matrices  A and B,  and on  entry  with
# >           TRANS = 'C' or 'c',  K  specifies  the number of rows of the
# >           matrices  A and B.  K must be at least zero.
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
# > \param[in] B
# > \verbatim
# >          B is COMPLEX array, dimension ( LDB, kb ), where kb is
# >           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
# >           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
# >           part of the array  B  must contain the matrix  B,  otherwise
# >           the leading  k by n  part of the array  B  must contain  the
# >           matrix B.
# > \endverbatim
# >
# > \param[in] LDB
# > \verbatim
# >          LDB is INTEGER
# >           On entry, LDB specifies the first dimension of B as declared
# >           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
# >           then  LDB must be at least  max( 1, n ), otherwise  LDB must
# >           be at least  max( 1, k ).
# > \endverbatim
# >
# > \param[in] BETA
# > \verbatim
# >          BETA is REAL
# >           On entry, BETA specifies the scalar beta.
# > \endverbatim
# >
# > \param[in,out] C
# > \verbatim
# >          C is COMPLEX array, dimension ( LDC, N )
# >           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
# >           upper triangular part of the array C must contain the upper
# >           triangular part  of the  hermitian matrix  and the strictly
# >           lower triangular part of C is not referenced.  On exit, the
# >           upper triangular part of the array  C is overwritten by the
# >           upper triangular part of the updated matrix.
# >           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
# >           lower triangular part of the array C must contain the lower
# >           triangular part  of the  hermitian matrix  and the strictly
# >           upper triangular part of C is not referenced.  On exit, the
# >           lower triangular part of the array  C is overwritten by the
# >           lower triangular part of the updated matrix.
# >           Note that the imaginary parts of the diagonal elements need
# >           not be set,  they are assumed to be zero,  and on exit they
# >           are set to zero.
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
# >
# >  -- Modified 8-Nov-93 to set C[J,J] to REAL( C[J,J] ) when BETA = 1.
# >     Ed Anderson, Cray Research Inc.
# > \endverbatim
# >
#  =====================================================================
from util import lsame
from xerbla import xerbla


def cher2k(UPLO, TRANS, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC):
    #
    #  -- Reference BLAS level3 routine (version 3.7.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     December 2016
    #
    #     .. Scalar Arguments ..
    #   COMPLEX ALPHA
    #   REAL BETA
    #   INTEGER K,LDA,LDB,LDC,N
    #   CHARACTER TRANS,UPLO
    #     ..
    #     .. Array Arguments ..
    #   COMPLEX A(LDA,*),B(LDB,*),C(LDC,*)
    #     ..
    #
    #  =====================================================================

    # Test the input parameters.
    if lsame(TRANS, "N"):
        NROWA = N
    else:
        NROWA = K
    UPPER = lsame(UPLO, "U")

    INFO = 0
    if (not UPPER) and (not lsame(UPLO, "L")):
        INFO = 1
    elif (not lsame(TRANS, "N")) and (not lsame(TRANS, "C")):
        INFO = 2
    elif N < 0:
        INFO = 3
    elif K < 0:
        INFO = 4
    elif LDA < max(1, NROWA):
        INFO = 7
    elif LDB < max(1, NROWA):
        INFO = 9
    elif LDC < max(1, N):
        INFO = 12
    if INFO != 0:
        xerbla("CHER2K", INFO)
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
                    for I in range(J - 1):
                        C[I, J] = BETA * C[I, J]
                    C[J, J] = BETA * C[J, J].real
        else:
            if BETA == 0:
                for J in range(N):
                    for I in range(J - 1, N):
                        C[I, J] = 0
            else:
                for J in range(N):
                    C[J, J] = BETA * C[J, J].real
                    for I in range(J, N):
                        C[I, J] = BETA * C[I, J]
        return

    # Start the operations.
    if lsame(TRANS, "N"):
        #
        #        Form  C := alpha*A*B**H + conjg( alpha )*B*A**H +
        #                   C.
        #
        if UPPER:
            for J in range(N):
                if BETA == 0:
                    for I in range(J):
                        C[I, J] = 0
                elif BETA != 1:
                    for I in range(J - 1):
                        C[I, J] = BETA * C[I, J]
                    C[J, J] = BETA * C[J, J].real
                else:
                    C[J, J] = C[J, J].real
                for L in range(K):
                    if (A[J, L] != 0) or (B[J, L] != 0):
                        TEMP1 = ALPHA * B[J, L].conjugate()
                        TEMP2 = (ALPHA * A[J, L]).conjugate()
                        for I in range(J - 1):
                            C[I, J] = C[I, J] + A[I, L] * TEMP1 + B[I, L] * TEMP2
                        C[J, J] = (
                            C[J, J].real + (A[J, L] * TEMP1 + B[J, L] * TEMP2).real
                        )
        else:
            for J in range(N):
                if BETA == 0:
                    for I in range(J - 1, N):
                        C[I, J] = 0
                elif BETA != 1:
                    for I in range(J, N):
                        C[I, J] = BETA * C[I, J]
                    C[J, J] = BETA * C[J, J].real
                else:
                    C[J, J] = C[J, J].real
                for L in range(K):
                    if (A[J, L] != 0) or (B[J, L] != 0):
                        TEMP1 = ALPHA * B[J, L].conjugate()
                        TEMP2 = (ALPHA * A[J, L]).conjugate()
                        for I in range(J, N):
                            C[I, J] = C[I, J] + A[I, L] * TEMP1 + B[I, L] * TEMP2
                        C[J, J] = (
                            C[J, J].real + (A[J, L] * TEMP1 + B[J, L] * TEMP2).real
                        )
    else:
        #
        #        Form  C := alpha*A**H*B + conjg( alpha )*B**H*A +
        #                   C.
        #
        if UPPER:
            for J in range(N):
                for I in range(J):
                    TEMP1 = 0
                    TEMP2 = 0
                    for L in range(K):
                        TEMP1 += A[L, I].conjugate() * B[L, J]
                        TEMP2 += (B(L, I)).conjugate() * A[L, J]
                    if I == J:
                        if BETA == 0:
                            C[J, J] = (ALPHA * TEMP1 + (ALPHA).conjugate() * TEMP2).real
                        else:
                            C[J, J] = (
                                BETA * C[J, J].real
                                + (ALPHA * TEMP1 + (ALPHA).conjugate() * TEMP2).real
                            )
                    else:
                        if BETA == 0:
                            C[I, J] = ALPHA * TEMP1 + (ALPHA).conjugate() * TEMP2
                        else:
                            C[I, J] = (
                                BETA * C[I, J]
                                + ALPHA * TEMP1
                                + (ALPHA).conjugate() * TEMP2
                            )
        else:
            for J in range(N):
                for I in range(J - 1, N):
                    TEMP1 = 0
                    TEMP2 = 0
                    for L in range(K):
                        TEMP1 += A[L, I].conjugate() * B[L, J]
                        TEMP2 += (B(L, I)).conjugate() * A[L, J]
                    if I == J:
                        if BETA == 0:
                            C[J, J] = (ALPHA * TEMP1 + (ALPHA).conjugate() * TEMP2).real
                        else:
                            C[J, J] = (
                                BETA * C[J, J].real
                                + (ALPHA * TEMP1 + (ALPHA).conjugate() * TEMP2).real
                            )
                    else:
                        if BETA == 0:
                            C[I, J] = ALPHA * TEMP1 + (ALPHA).conjugate() * TEMP2
                        else:
                            C[I, J] = (
                                BETA * C[I, J]
                                + ALPHA * TEMP1
                                + (ALPHA).conjugate() * TEMP2
                            )
