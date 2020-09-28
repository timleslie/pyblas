# > \brief \b CHEMM
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def CHEMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
#
#       .. Scalar Arguments ..
#       COMPLEX ALPHA,BETA
#       INTEGER LDA,LDB,LDC,M,N
#       CHARACTER SIDE,UPLO
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
# > CHEMM  performs one of the matrix-matrix operations
# >
# >    C := alpha*A*B + beta*C,
# >
# > or
# >
# >    C := alpha*B*A + beta*C,
# >
# > where alpha and beta are scalars, A is an hermitian matrix and  B and
# > C are m by n matrices.
# > \endverbatim
#
#  Arguments:
#  ==========
#
# > \param[in] SIDE
# > \verbatim
# >          SIDE is CHARACTER*1
# >           On entry,  SIDE  specifies whether  the  hermitian matrix  A
# >           appears on the  left or right  in the  operation as follows:
# >
# >              SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,
# >
# >              SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,
# > \endverbatim
# >
# > \param[in] UPLO
# > \verbatim
# >          UPLO is CHARACTER*1
# >           On  entry,   UPLO  specifies  whether  the  upper  or  lower
# >           triangular  part  of  the  hermitian  matrix   A  is  to  be
# >           referenced as follows:
# >
# >              UPLO = 'U' or 'u'   Only the upper triangular part of the
# >                                  hermitian matrix is to be referenced.
# >
# >              UPLO = 'L' or 'l'   Only the lower triangular part of the
# >                                  hermitian matrix is to be referenced.
# > \endverbatim
# >
# > \param[in] M
# > \verbatim
# >          M is INTEGER
# >           On entry,  M  specifies the number of rows of the matrix  C.
# >           M  must be at least zero.
# > \endverbatim
# >
# > \param[in] N
# > \verbatim
# >          N is INTEGER
# >           On entry, N specifies the number of columns of the matrix C.
# >           N  must be at least zero.
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
# >           m  when  SIDE = 'L' or 'l'  and is n  otherwise.
# >           Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of
# >           the array  A  must contain the  hermitian matrix,  such that
# >           when  UPLO = 'U' or 'u', the leading m by m upper triangular
# >           part of the array  A  must contain the upper triangular part
# >           of the  hermitian matrix and the  strictly  lower triangular
# >           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
# >           the leading  m by m  lower triangular part  of the  array  A
# >           must  contain  the  lower triangular part  of the  hermitian
# >           matrix and the  strictly upper triangular part of  A  is not
# >           referenced.
# >           Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of
# >           the array  A  must contain the  hermitian matrix,  such that
# >           when  UPLO = 'U' or 'u', the leading n by n upper triangular
# >           part of the array  A  must contain the upper triangular part
# >           of the  hermitian matrix and the  strictly  lower triangular
# >           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
# >           the leading  n by n  lower triangular part  of the  array  A
# >           must  contain  the  lower triangular part  of the  hermitian
# >           matrix and the  strictly upper triangular part of  A  is not
# >           referenced.
# >           Note that the imaginary parts  of the diagonal elements need
# >           not be set, they are assumed to be zero.
# > \endverbatim
# >
# > \param[in] LDA
# > \verbatim
# >          LDA is INTEGER
# >           On entry, LDA specifies the first dimension of A as declared
# >           in the  calling (sub) program. When  SIDE = 'L' or 'l'  then
# >           LDA must be at least  max( 1, m ), otherwise  LDA must be at
# >           least max( 1, n ).
# > \endverbatim
# >
# > \param[in] B
# > \verbatim
# >          B is COMPLEX array, dimension ( LDB, N )
# >           Before entry, the leading  m by n part of the array  B  must
# >           contain the matrix B.
# > \endverbatim
# >
# > \param[in] LDB
# > \verbatim
# >          LDB is INTEGER
# >           On entry, LDB specifies the first dimension of B as declared
# >           in  the  calling  (sub)  program.   LDB  must  be  at  least
# >           max( 1, m ).
# > \endverbatim
# >
# > \param[in] BETA
# > \verbatim
# >          BETA is COMPLEX
# >           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
# >           supplied as zero then C need not be set on input.
# > \endverbatim
# >
# > \param[in,out] C
# > \verbatim
# >          C is COMPLEX array, dimension ( LDC, N )
# >           Before entry, the leading  m by n  part of the array  C must
# >           contain the matrix  C,  except when  beta  is zero, in which
# >           case C need not be set on entry.
# >           On exit, the array  C  is overwritten by the  m by n updated
# >           matrix.
# > \endverbatim
# >
# > \param[in] LDC
# > \verbatim
# >          LDC is INTEGER
# >           On entry, LDC specifies the first dimension of C as declared
# >           in  the  calling  (sub)  program.   LDC  must  be  at  least
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


def CHEMM(SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB, BETA, C, LDC):
    #
    #  -- Reference BLAS level3 routine (version 3.7.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     December 2016
    #
    #     .. Scalar Arguments ..
    #   COMPLEX ALPHA,BETA
    #   INTEGER LDA,LDB,LDC,M,N
    #   CHARACTER SIDE,UPLO
    #     ..
    #     .. Array Arguments ..
    #   COMPLEX A(LDA,*),B(LDB,*),C(LDC,*)
    #     ..
    #
    #  =====================================================================
    #
    #     Set NROWA as the number of rows of A.
    #
    if lsame(SIDE, "L"):
        NROWA = M
    else:
        NROWA = N
    UPPER = lsame(UPLO, "U")

    # Test the input parameters.
    INFO = 0
    if (not lsame(SIDE, "L")) and (not lsame(SIDE, "R")):
        INFO = 1
    elif (not UPPER) and (not lsame(UPLO, "L")):
        INFO = 2
    elif M < 0:
        INFO = 3
    elif N < 0:
        INFO = 4
    elif LDA < max(1, NROWA):
        INFO = 7
    elif LDB < max(1, M):
        INFO = 9
    elif LDC < max(1, M):
        INFO = 12
    if INFO != 0:
        xerbla("CHEMM ", INFO)
        return

    # Quick return if possible.
    if (M == 0) or (N == 0) or ((ALPHA == 0) and (BETA == 1)):
        return

    # And when  alpha==zero.
    if ALPHA == 0:
        if BETA == 0:
            for J in range(N):
                for I in range(M):
                    C[I, J] = 0
        else:
            for J in range(N):
                for I in range(M):
                    C[I, J] = BETA * C[I, J]
        return

    # Start the operations.
    if lsame(SIDE, "L"):
        #
        #        Form  C := alpha*A*B + beta*C.
        #
        if UPPER:
            for J in range(N):
                for I in range(M):
                    TEMP1 = ALPHA * B[I, J]
                    TEMP2 = 0
                    for K in range(I - 1):
                        C[K, J] = C[K, J] + TEMP1 * A[K, I]
                        TEMP2 += B[K, J] * A[K, I].conjugate()
                    if BETA == 0:
                        C[I, J] = TEMP1 * A[I, I].real + ALPHA * TEMP2
                    else:
                        C[I, J] = BETA * C[I, J] + TEMP1 * A[I, I].real + ALPHA * TEMP2
        else:
            for J in range(N):
                for I in range(M - 1, -1, -1):
                    TEMP1 = ALPHA * B[I, J]
                    TEMP2 = 0
                    for K in range(I, M):
                        C[K, J] = C[K, J] + TEMP1 * A[K, I]
                        TEMP2 += B[K, J] * A[K, I].conjugate()
                    if BETA == 0:
                        C[I, J] = TEMP1 * A[I, I].real + ALPHA * TEMP2
                    else:
                        C[I, J] = BETA * C[I, J] + TEMP1 * A[I, I].real + ALPHA * TEMP2
    else:
        #
        #        Form  C := alpha*B*A + beta*C.
        #
        for J in range(N):
            TEMP1 = ALPHA * A[J, J].real
            if BETA == 0:
                for I in range(M):
                    C[I, J] = TEMP1 * B[I, J]
            else:
                for I in range(M):
                    C[I, J] = BETA * C[I, J] + TEMP1 * B[I, J]
            for K in range(J - 1):
                if UPPER:
                    TEMP1 = ALPHA * A[K, J]
                else:
                    TEMP1 = ALPHA * A[J, K].conjugate()
                for I in range(M):
                    C[I, J] = C[I, J] + TEMP1 * B[I, K]
            for K in range(J, N):
                if UPPER:
                    TEMP1 = ALPHA * A[J, K].conjugate()
                else:
                    TEMP1 = ALPHA * A[K, J]
                for I in range(M):
                    C[I, J] = C[I, J] + TEMP1 * B[I, K]
