# > \brief \b CHPR2
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def CHPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP)
#
#       .. Scalar Arguments ..
#       COMPLEX ALPHA
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
# > CHPR2  performs the hermitian rank 2 operation
# >
# >    A := alpha*x*y**H + conjg( alpha )*y*x**H + A,
# >
# > where alpha is a scalar, x and y are n element vectors and A is an
# > n by n hermitian matrix, supplied in packed form.
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
# > \param[in] Y
# > \verbatim
# >          Y is COMPLEX array, dimension at least
# >           ( 1 + ( n - 1 )*abs( INCY ) ).
# >           Before entry, the incremented array Y must contain the n
# >           element vector y.
# > \endverbatim
# >
# > \param[in] INCY
# > \verbatim
# >          INCY is INTEGER
# >           On entry, INCY specifies the increment for the elements of
# >           Y. INCY must not be zero.
# > \endverbatim
# >
# > \param[in,out] AP
# > \verbatim
# >          AP is COMPLEX array, dimension at least
# >           ( ( n*( n + 1 ) )/2 ).
# >           Before entry with  UPLO = 'U' or 'u', the array AP must
# >           contain the upper triangular part of the hermitian matrix
# >           packed sequentially, column by column, so that AP( 1 )
# >           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
# >           and a( 2, 2 ) respectively, and so on. On exit, the array
# >           AP is overwritten by the upper triangular part of the
# >           updated matrix.
# >           Before entry with UPLO = 'L' or 'l', the array AP must
# >           contain the lower triangular part of the hermitian matrix
# >           packed sequentially, column by column, so that AP( 1 )
# >           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
# >           and a( 3, 1 ) respectively, and so on. On exit, the array
# >           AP is overwritten by the lower triangular part of the
# >           updated matrix.
# >           Note that the imaginary parts of the diagonal elements need
# >           not be set, they are assumed to be zero, and on exit they
# >           are set to zero.
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


def chpr2(UPLO, N, ALPHA, X, INCX, Y, INCY, AP):
    #
    #  -- Reference BLAS level2 routine (version 3.7.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     December 2016
    #
    #     .. Scalar Arguments ..
    #   COMPLEX ALPHA
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
        INFO = 5
    elif INCY == 0:
        INFO = 7
    if INFO != 0:
        xerbla("CHPR2 ", INFO)

    # Quick return if possible.
    if (N == 0) or (ALPHA == 0):
        return

    # Set up the start points in X and Y if the increments are not both unity.
    if (INCX != 1) or (INCY != 1):
        if INCX > 0:
            KX = 1
        else:
            KX = 1 - (N - 1) * INCX
        if INCY > 0:
            KY = 1
        else:
            KY = 1 - (N - 1) * INCY
        JX = KX
        JY = KY

    # Start the operations. In this version the elements of the array AP
    # are accessed sequentially with one pass through AP.
    KK = 1
    if lsame(UPLO, "U"):
        # Form  A  when upper triangle is stored in AP.
        if (INCX == 1) and (INCY == 1):
            for J in range(N):
                if (X[J] != 0) or (Y[J] != 0):
                    TEMP1 = ALPHA * (Y[J]).conjugate()
                    TEMP2 = (ALPHA * X[J]).conjugate()
                    K = KK
                    for I in range(J - 1):
                        AP[K] = AP[K] + X[I] * TEMP1 + Y[I] * TEMP2
                        K += 1
                    AP[KK + J - 1] = (AP[KK + J - 1]).real + (
                        X[J] * TEMP1 + Y[J] * TEMP2
                    ).real
                else:
                    AP[KK + J - 1] = (AP[KK + J - 1]).real
                KK += J
        else:
            for J in range(N):
                if (X[JX] != 0) or (Y[JY] != 0):
                    TEMP1 = ALPHA * (Y[JY]).conjugate()
                    TEMP2 = (ALPHA * X[JX]).conjugate()
                    IX = KX
                    IY = KY
                    for K in range(KK - 1, KK + J - 2):
                        AP[K] = AP[K] + X[IX] * TEMP1 + Y[IY] * TEMP2
                        IX += INCX
                        IY += INCY
                    AP[KK + J - 1] = (AP[KK + J - 1]).real + (
                        X[JX] * TEMP1 + Y[JY] * TEMP2
                    ).real
                else:
                    AP[KK + J - 1] = (AP[KK + J - 1]).real
                JX += INCX
                JY += INCY
                KK += J
    else:
        # Form  A  when lower triangle is stored in AP.
        if (INCX == 1) and (INCY == 1):
            for J in range(N):
                if (X[J] != 0) or (Y[J] != 0):
                    TEMP1 = ALPHA * (Y[J]).conjugate()
                    TEMP2 = (ALPHA * X[J]).conjugate()
                    AP[KK] = (AP[KK]).real + (X[J] * TEMP1 + Y[J] * TEMP2).real
                    K = KK + 1
                    for I in range(J, N):
                        AP[K] = AP[K] + X[I] * TEMP1 + Y[I] * TEMP2
                        K += 1
                else:
                    AP[KK] = (AP[KK]).real
                KK = KK + N - J + 1
        else:
            for J in range(N):
                if (X[JX] != 0) or (Y[JY] != 0):
                    TEMP1 = ALPHA * (Y[JY]).conjugate()
                    TEMP2 = (ALPHA * X[JX]).conjugate()
                    AP[KK] = (AP[KK]).real + (X[JX] * TEMP1 + Y[JY] * TEMP2).real
                    IX = JX
                    IY = JY
                    for K in range(KK, KK + N - J):
                        IX += INCX
                        IY += INCY
                        AP[K] = AP[K] + X[IX] * TEMP1 + Y[IY] * TEMP2
                else:
                    AP[KK] = (AP[KK]).real
                JX += INCX
                JY += INCY
                KK = KK + N - J + 1
