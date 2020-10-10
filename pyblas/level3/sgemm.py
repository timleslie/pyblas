# > \brief \b SGEMM
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def SGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
#
#       .. Scalar Arguments ..
#       REAL ALPHA,BETA
#       INTEGER K,LDA,LDB,LDC,M,N
#       CHARACTER TRANSA,TRANSB
#       ..
#       .. Array Arguments ..
#       REAL A(LDA,*),B(LDB,*),C(LDC,*)
#       ..
#
#
# > \par Purpose:
#  =============
# >
# > \verbatim
# >
# > SGEMM  performs one of the matrix-matrix operations
# >
# >    C := alpha*op( A )*op( B ) + beta*C,
# >
# > where  op( X ) is one of
# >
# >    op( X ) = X   or   op( X ) = X**T,
# >
# > alpha and beta are scalars, and A, B and C are matrices, with op( A )
# > an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
# > \endverbatim
#
#  Arguments:
#  ==========
#
# > \param[in] TRANSA
# > \verbatim
# >          TRANSA is CHARACTER*1
# >           On entry, TRANSA specifies the form of op( A ) to be used in
# >           the matrix multiplication as follows:
# >
# >              TRANSA = 'N' or 'n',  op( A ) = A.
# >
# >              TRANSA = 'T' or 't',  op( A ) = A**T.
# >
# >              TRANSA = 'C' or 'c',  op( A ) = A**T.
# > \endverbatim
# >
# > \param[in] TRANSB
# > \verbatim
# >          TRANSB is CHARACTER*1
# >           On entry, TRANSB specifies the form of op( B ) to be used in
# >           the matrix multiplication as follows:
# >
# >              TRANSB = 'N' or 'n',  op( B ) = B.
# >
# >              TRANSB = 'T' or 't',  op( B ) = B**T.
# >
# >              TRANSB = 'C' or 'c',  op( B ) = B**T.
# > \endverbatim
# >
# > \param[in] M
# > \verbatim
# >          M is INTEGER
# >           On entry,  M  specifies  the number  of rows  of the  matrix
# >           op( A )  and of the  matrix  C.  M  must  be at least  zero.
# > \endverbatim
# >
# > \param[in] N
# > \verbatim
# >          N is INTEGER
# >           On entry,  N  specifies the number  of columns of the matrix
# >           op( B ) and the number of columns of the matrix C. N must be
# >           at least zero.
# > \endverbatim
# >
# > \param[in] K
# > \verbatim
# >          K is INTEGER
# >           On entry,  K  specifies  the number of columns of the matrix
# >           op( A ) and the number of rows of the matrix op( B ). K must
# >           be at least  zero.
# > \endverbatim
# >
# > \param[in] ALPHA
# > \verbatim
# >          ALPHA is REAL
# >           On entry, ALPHA specifies the scalar alpha.
# > \endverbatim
# >
# > \param[in] A
# > \verbatim
# >          A is REAL array, dimension ( LDA, ka ), where ka is
# >           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
# >           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
# >           part of the array  A  must contain the matrix  A,  otherwise
# >           the leading  k by m  part of the array  A  must contain  the
# >           matrix A.
# > \endverbatim
# >
# > \param[in] LDA
# > \verbatim
# >          LDA is INTEGER
# >           On entry, LDA specifies the first dimension of A as declared
# >           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
# >           LDA must be at least  max( 1, m ), otherwise  LDA must be at
# >           least  max( 1, k ).
# > \endverbatim
# >
# > \param[in] B
# > \verbatim
# >          B is REAL array, dimension ( LDB, kb ), where kb is
# >           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
# >           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
# >           part of the array  B  must contain the matrix  B,  otherwise
# >           the leading  n by k  part of the array  B  must contain  the
# >           matrix B.
# > \endverbatim
# >
# > \param[in] LDB
# > \verbatim
# >          LDB is INTEGER
# >           On entry, LDB specifies the first dimension of B as declared
# >           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
# >           LDB must be at least  max( 1, k ), otherwise  LDB must be at
# >           least  max( 1, n ).
# > \endverbatim
# >
# > \param[in] BETA
# > \verbatim
# >          BETA is REAL
# >           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
# >           supplied as zero then C need not be set on input.
# > \endverbatim
# >
# > \param[in,out] C
# > \verbatim
# >          C is REAL array, dimension ( LDC, N )
# >           Before entry, the leading  m by n  part of the array  C must
# >           contain the matrix  C,  except when  beta  is zero, in which
# >           case C need not be set on entry.
# >           On exit, the array  C  is overwritten by the  m by n  matrix
# >           ( alpha*op( A )*op( B ) + beta*C ).
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
# > \ingroup single_blas_level3
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


def SGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC):
    #
    #  -- Reference BLAS level3 routine (version 3.7.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     December 2016
    #
    #     .. Scalar Arguments ..
    #   REAL ALPHA,BETA
    #   INTEGER K,LDA,LDB,LDC,M,N
    #   CHARACTER TRANSA,TRANSB
    #     ..
    #     .. Array Arguments ..
    #   REAL A(LDA,*),B(LDB,*),C(LDC,*)
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
    #   REAL TEMP
    #   INTEGER I,INFO,J,L,NCOLA,NROWA,NROWB
    #   LOGICAL NOTA,NOTB
    #     ..
    #     .. Parameters ..
    #   REAL ONE,ZERO
    #   PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
    #     ..
    #
    #     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
    #     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
    #     and  columns of  A  and the  number of  rows  of  B  respectively.
    #
    NOTA = lsame(TRANSA, "N")
    NOTB = lsame(TRANSB, "N")
    if NOTA:
        NROWA = M
    else:
        NROWA = K
    if NOTB:
        NROWB = K
    else:
        NROWB = N

    # Test the input parameters.
    INFO = 0
    if (not NOTA) and (not lsame(TRANSA, "C")) and (not lsame(TRANSA, "T")):
        INFO = 1
    elif (not NOTB) and (not lsame(TRANSB, "C")) and (not lsame(TRANSB, "T")):
        INFO = 2
    elif M < 0:
        INFO = 3
    elif N < 0:
        INFO = 4
    elif K < 0:
        INFO = 5
    elif LDA < max(1, NROWA):
        INFO = 8
    elif LDB < max(1, NROWB):
        INFO = 10
    elif LDC < max(1, M):
        INFO = 13
    if INFO != 0:
        xerbla("SGEMM ", INFO)
        return

    # Quick return if possible.
    if (M == 0) or (N == 0) or (((ALPHA == 0) or (K == 0)) and (BETA == 1)):
        return
    #
    #     And if  alpha==zero.
    #
    if ALPHA == 0:
        if BETA == 0:
            for J in range(N):
                for I in range(M):
                    C[I, J] = 0
        else:
            for J in range(N):
                for I in range(M):
                    C[I, J] *= BETA
        return

    # Start the operations.
    if NOTB:
        if NOTA:
            #
            #           Form  C := alpha*A*B + beta*C.
            #
            for J in range(N):
                if BETA == 0:
                    for I in range(M):
                        C[I, J] = 0
                elif BETA != 1:
                    for I in range(M):
                        C[I, J] *= BETA
                for L in range(K):
                    TEMP = ALPHA * B[L, J]
                    for I in range(M):
                        C[I, J] += TEMP * A[I, L]
        else:
            #
            #           Form  C := alpha*A**T*B + beta*C
            #
            for J in range(N):
                for I in range(M):
                    TEMP = 0
                    for L in range(K):
                        TEMP += A[L, I] * B[L, J]
                    if BETA == 0:
                        C[I, J] = ALPHA * TEMP
                    else:
                        C[I, J] = ALPHA * TEMP + BETA * C[I, J]
    else:
        if NOTA:
            #
            #           Form  C := alpha*A*B**T + beta*C
            #
            for J in range(N):
                if BETA == 0:
                    for I in range(M):
                        C[I, J] = 0
                elif BETA != 1:
                    for I in range(M):
                        C[I, J] *= BETA
                for L in range(K):
                    TEMP = ALPHA * B[J, L]
                    for I in range(M):
                        C[I, J] += TEMP * A[I, L]
        else:
            # Form  C := alpha*A**T*B**T + beta*C
            for J in range(N):
                for I in range(M):
                    TEMP = 0
                    for L in range(K):
                        TEMP += A[L, I] * B[J, L]
                    if BETA == 0:
                        C[I, J] = ALPHA * TEMP
                    else:
                        C[I, J] = ALPHA * TEMP + BETA * C[I, J]
