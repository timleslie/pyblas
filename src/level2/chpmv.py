# > \brief \b CHPMV
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def CHPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)
#
#       .. Scalar Arguments ..
#       COMPLEX ALPHA,BETA
#       INTEGER INCX,INCY,N
#       CHARACTER UPLO
#       ..
#       .. Array Arguments ..
#       COMPLEX AP(*),X(*),Y(*)
#       ..
#
#
# > \par Purpose:
#  =============
# >
# > \verbatim
# >
# > CHPMV  performs the matrix-vector operation
# >
# >    y := alpha*A*x + beta*y,
# >
# > where alpha and beta are scalars, x and y are n element vectors and
# > A is an n by n hermitian matrix, supplied in packed form.
# > \endverbatim
#
#  Arguments:
#  ==========
#
# > \param[in] UPLO
# > \verbatim
# >          UPLO is CHARACTER*1
# >           On entry, UPLO specifies whether the upper or lower
# >           triangular part of the matrix A is supplied in the packed
# >           array AP as follows:
# >
# >              UPLO = 'U' or 'u'   The upper triangular part of A is
# >                                  supplied in AP.
# >
# >              UPLO = 'L' or 'l'   The lower triangular part of A is
# >                                  supplied in AP.
# > \endverbatim
# >
# > \param[in] N
# > \verbatim
# >          N is INTEGER
# >           On entry, N specifies the order of the matrix A.
# >           N must be at least zero.
# > \endverbatim
# >
# > \param[in] ALPHA
# > \verbatim
# >          ALPHA is COMPLEX
# >           On entry, ALPHA specifies the scalar alpha.
# > \endverbatim
# >
# > \param[in] AP
# > \verbatim
# >          AP is COMPLEX array, dimension at least
# >           ( ( n*( n + 1 ) )/2 ).
# >           Before entry with UPLO = 'U' or 'u', the array AP must
# >           contain the upper triangular part of the hermitian matrix
# >           packed sequentially, column by column, so that AP( 1 )
# >           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
# >           and a( 2, 2 ) respectively, and so on.
# >           Before entry with UPLO = 'L' or 'l', the array AP must
# >           contain the lower triangular part of the hermitian matrix
# >           packed sequentially, column by column, so that AP( 1 )
# >           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
# >           and a( 3, 1 ) respectively, and so on.
# >           Note that the imaginary parts of the diagonal elements need
# >           not be set and are assumed to be zero.
# > \endverbatim
# >
# > \param[in] X
# > \verbatim
# >          X is COMPLEX array, dimension at least
# >           ( 1 + ( n - 1 )*abs( INCX ) ).
# >           Before entry, the incremented array X must contain the n
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
# >           ( 1 + ( n - 1 )*abs( INCY ) ).
# >           Before entry, the incremented array Y must contain the n
# >           element vector y. On exit, Y is overwritten by the updated
# >           vector y.
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


def chpmv(UPLO, N, ALPHA, AP, X, INCX, BETA, Y, INCY):
    #
    #  -- Reference BLAS level2 routine (version 3.7.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     December 2016
    #
    #     .. Scalar Arguments ..
    #   COMPLEX ALPHA,BETA
    #   INTEGER INCX,INCY,N
    #   CHARACTER UPLO
    #     ..
    #     .. Array Arguments ..
    #   COMPLEX AP(*),X(*),Y(*)
    #     ..
    #
    #  =====================================================================

    # Test the input parameters.
    INFO = 0
    if not lsame(UPLO, "U") and not lsame(UPLO, "L"):
        INFO = 1
    elif N < 0:
        INFO = 2
    elif INCX == 0:
        INFO = 6
    elif INCY == 0:
        INFO = 9
    if INFO != 0:
        xerbla("CHPMV ", INFO)

    # Quick return if possible.
    if (N == 0) or ((ALPHA == 0) and (BETA == 1)):
        return

    # Set up the start points in  X  and  Y.
    if INCX > 0:
        KX = 1
    else:
        KX = 1 - (N - 1) * INCX
    if INCY > 0:
        KY = 1
    else:
        KY = 1 - (N - 1) * INCY

    # Start the operations. In this version the elements of the array AP
    # are accessed sequentially with one pass through AP.

    # First form  y := beta*y.
    if BETA != 1:
        if INCY == 1:
            if BETA == 0:
                for I in range(N):
                    Y[I] = 0
            else:
                for I in range(N):
                    Y[I] = BETA * Y[I]
        else:
            IY = KY
            if BETA == 0:
                for I in range(N):
                    Y[IY] = 0
                    IY += INCY
            else:
                for I in range(N):
                    Y[IY] = BETA * Y[IY]
                    IY += INCY
    if ALPHA == 0:
        return
    KK = 1
    if lsame(UPLO, "U"):
        # Form  y  when AP contains the upper triangle.
        if (INCX == 1) and (INCY == 1):
            for J in range(N):
                TEMP1 = ALPHA * X[J]
                TEMP2 = 0
                K = KK
                for I in range(J - 1):
                    Y[I] = Y[I] + TEMP1 * AP[K]
                    TEMP2 += (AP[K]).conjugate() * X[I]
                    K += 1
                Y[J] = Y[J] + TEMP1 * (AP[KK + J - 1]).real + ALPHA * TEMP2
                KK += J
        else:
            JX = KX
            JY = KY
            for J in range(N):
                TEMP1 = ALPHA * X[JX]
                TEMP2 = 0
                IX = KX
                IY = KY
                for K in range(KK - 1, KK + J - 2):
                    Y[IY] = Y[IY] + TEMP1 * AP[K]
                    TEMP2 += (AP[K]).conjugate() * X[IX]
                    IX += INCX
                    IY += INCY
                Y[JY] += TEMP1 * (AP[KK + J - 1]).real + ALPHA * TEMP2
                JX += INCX
                JY += INCY
                KK += J
    else:
        # Form  y  when AP contains the lower triangle.
        if (INCX == 1) and (INCY == 1):
            for J in range(N):
                TEMP1 = ALPHA * X[J]
                TEMP2 = 0
                Y[J] = Y[J] + TEMP1 * (AP[KK]).real
                K = KK + 1
                for I in range(J, N):
                    Y[I] = Y[I] + TEMP1 * AP[K]
                    TEMP2 += (AP[K]).conjugate() * X[I]
                    K += 1
                Y[J] = Y[J] + ALPHA * TEMP2
                KK = KK + (N - J + 1)
        else:
            JX = KX
            JY = KY
            for J in range(N):
                TEMP1 = ALPHA * X[JX]
                TEMP2 = 0
                Y[JY] += TEMP1 * (AP[KK]).real
                IX = JX
                IY = JY
                for K in range(KK, KK + N - J):
                    IX += INCX
                    IY += INCY
                    Y[IY] = Y[IY] + TEMP1 * AP[K]
                    TEMP2 += (AP[K]).conjugate() * X[IX]
                Y[JY] += ALPHA * TEMP2
                JX += INCX
                JY += INCY
                KK = KK + (N - J + 1)
