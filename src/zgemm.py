# > \brief \b ZGEMM
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
#
#       .. Scalar Arguments ..
#       COMPLEX*16 ALPHA,BETA
#       INTEGER K,LDA,LDB,LDC,M,N
#       CHARACTER TRANSA,TRANSB
#       ..
#       .. Array Arguments ..
#       COMPLEX*16 A(LDA,*),B(LDB,*),C(LDC,*)
#       ..
#
#
# > \par Purpose:
#  =============
# >
# > \verbatim
# >
# > ZGEMM  performs one of the matrix-matrix operations
# >
# >    C := alpha*op( A )*op( B ) + beta*C,
# >
# > where  op( X ) is one of
# >
# >    op( X ) = X   or   op( X ) = X**T   or   op( X ) = X**H,
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
# >              TRANSA = 'C' or 'c',  op( A ) = A**H.
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
# >              TRANSB = 'C' or 'c',  op( B ) = B**H.
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
# >          ALPHA is COMPLEX*16
# >           On entry, ALPHA specifies the scalar alpha.
# > \endverbatim
# >
# > \param[in] A
# > \verbatim
# >          A is COMPLEX*16 array, dimension ( LDA, ka ), where ka is
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
# >          B is COMPLEX*16 array, dimension ( LDB, kb ), where kb is
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
# >          BETA is COMPLEX*16
# >           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
# >           supplied as zero then C need not be set on input.
# > \endverbatim
# >
# > \param[in,out] C
# > \verbatim
# >          C is COMPLEX*16 array, dimension ( LDC, N )
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


def ZGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC):
    #
    #  -- Reference BLAS level3 routine (version 3.7.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     December 2016
    #
    #     .. Scalar Arguments ..
    #   COMPLEX*16 ALPHA,BETA
    #   INTEGER K,LDA,LDB,LDC,M,N
    #   CHARACTER TRANSA,TRANSB
    #     ..
    #     .. Array Arguments ..
    #   COMPLEX*16 A(LDA,*),B(LDB,*),C(LDC,*)
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
    #   INTRINSIC DCONJG,MAX
    #     ..
    #     .. Local Scalars ..
    #   COMPLEX*16 TEMP
    #   INTEGER I,INFO,J,L,NCOLA,NROWA,NROWB
    #   LOGICAL CONJA,CONJB,NOTA,NOTB
    #     ..
    #     .. Parameters ..
    #   COMPLEX*16 ONE
    #   PARAMETER (ONE= (1.0D+0,0.0D+0))
    #   COMPLEX*16 ZERO
    #   PARAMETER (ZERO= (0.0D+0,0.0D+0))
    #     ..
    #
    #     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
    #     conjugated or transposed, set  CONJA and CONJB  as true if  A  and
    #     B  respectively are to be  transposed but  not conjugated  and set
    #     NROWA, NCOLA and  NROWB  as the number of rows and  columns  of  A
    #     and the number of rows of  B  respectively.
    #
    NOTA = lsame(TRANSA, "N")
    NOTB = lsame(TRANSB, "N")
    CONJA = lsame(TRANSA, "C")
    CONJB = lsame(TRANSB, "C")
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
    if (not NOTA) and (not CONJA) and (not lsame(TRANSA, "T")):
        INFO = 1
    elif (not NOTB) and (not CONJB) and (not lsame(TRANSB, "T")):
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
        xerbla("ZGEMM ", INFO)

    # Quick return if possible.
    if (M == 0) or (N == 0) or (((ALPHA == 0) or (K == 0)) and (BETA == 1)):
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
                        C[I, J] = BETA * C[I, J]
                for L in range(K):
                    TEMP = ALPHA * B[L, J]
                    for I in range(M):
                        C[I, J] = C[I, J] + TEMP * A[I, L]
        elif CONJA:
            #
            #           Form  C := alpha*A**H*B + beta*C.
            #
            for J in range(N):
                for I in range(M):
                    TEMP = 0
                    for L in range(K):
                        TEMP += (A[L, I]).conjugate() * B[L, J]
                    if BETA == 0:
                        C[I, J] = ALPHA * TEMP
                    else:
                        C[I, J] = ALPHA * TEMP + BETA * C[I, J]
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
    elif NOTA:
        if CONJB:
            #
            #           Form  C := alpha*A*B**H + beta*C.
            #
            for J in range(N):
                if BETA == 0:
                    for I in range(M):
                        C[I, J] = 0
                elif BETA != 1:
                    for I in range(M):
                        C[I, J] = BETA * C[I, J]
                for L in range(K):
                    TEMP = ALPHA * (B[J, L]).conjugate()
                    for I in range(M):
                        C[I, J] = C[I, J] + TEMP * A[I, L]
        else:
            #
            #           Form  C := alpha*A*B**T + beta*C
            #
            for J in range(N):
                if BETA == 0:
                    for I in range(M):
                        C[I, J] = 0
                elif BETA != 1:
                    for I in range(M):
                        C[I, J] = BETA * C[I, J]
                for L in range(K):
                    TEMP = ALPHA * B[J, L]
                    for I in range(M):
                        C[I, J] = C[I, J] + TEMP * A[I, L]
    elif CONJA:
        if CONJB:
            #
            #           Form  C := alpha*A**H*B**H + beta*C.
            #
            for J in range(N):
                for I in range(M):
                    TEMP = 0
                    for L in range(K):
                        TEMP += (A[L, I]).conjugate() * (B[J, L]).conjugate()
                    if BETA == 0:
                        C[I, J] = ALPHA * TEMP
                    else:
                        C[I, J] = ALPHA * TEMP + BETA * C[I, J]
        else:
            #
            #           Form  C := alpha*A**H*B**T + beta*C
            #
            for J in range(N):
                for I in range(M):
                    TEMP = 0
                    for L in range(K):
                        TEMP += (A[L, I]).conjugate() * B[J, L]
                    if BETA == 0:
                        C[I, J] = ALPHA * TEMP
                    else:
                        C[I, J] = ALPHA * TEMP + BETA * C[I, J]
    else:
        if CONJB:
            #
            #           Form  C := alpha*A**T*B**H + beta*C
            #
            for J in range(N):
                for I in range(M):
                    TEMP = 0
                    for L in range(K):
                        TEMP += A[L, I] * (B[J, L]).conjugate()
                    if BETA == 0:
                        C[I, J] = ALPHA * TEMP
                    else:
                        C[I, J] = ALPHA * TEMP + BETA * C[I, J]
        else:
            #
            #           Form  C := alpha*A**T*B**T + beta*C
            #
            for J in range(N):
                for I in range(M):
                    TEMP = 0
                    for L in range(K):
                        TEMP += A[L, I] * B[J, L]
                    if BETA == 0:
                        C[I, J] = ALPHA * TEMP
                    else:
                        C[I, J] = ALPHA * TEMP + BETA * C[I, J]
