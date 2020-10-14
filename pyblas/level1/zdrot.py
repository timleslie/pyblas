# > \brief \b ZDROT
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def ZDROT( N, CX, INCX, CY, INCY, C, S )
#
#       .. Scalar Arguments ..
#       INTEGER            INCX, INCY, N
#       DOUBLE PRECISION   C, S
#       ..
#       .. Array Arguments ..
#       COMPLEX*16         CX( * ), CY( * )
#       ..
#
#
# > \par Purpose:
#  =============
# >
# > \verbatim
# >
# > Applies a plane rotation, where the cos and sin (c and s) are real
# > and the vectors cx and cy are complex.
# > jack dongarra, linpack, 3/11/78.
# > \endverbatim
#
#  Arguments:
#  ==========
#
# > \param[in] N
# > \verbatim
# >          N is INTEGER
# >           On entry, N specifies the order of the vectors cx and cy.
# >           N must be at least zero.
# > \endverbatim
# >
# > \param[in,out] CX
# > \verbatim
# >          CX is COMPLEX*16 array, dimension at least
# >           ( 1 + ( N - 1 )*abs( INCX ) ).
# >           Before entry, the incremented array CX must contain the n
# >           element vector cx. On exit, CX is overwritten by the updated
# >           vector cx.
# > \endverbatim
# >
# > \param[in] INCX
# > \verbatim
# >          INCX is INTEGER
# >           On entry, INCX specifies the increment for the elements of
# >           CX. INCX must not be zero.
# > \endverbatim
# >
# > \param[in,out] CY
# > \verbatim
# >          CY is COMPLEX*16 array, dimension at least
# >           ( 1 + ( N - 1 )*abs( INCY ) ).
# >           Before entry, the incremented array CY must contain the n
# >           element vector cy. On exit, CY is overwritten by the updated
# >           vector cy.
# > \endverbatim
# >
# > \param[in] INCY
# > \verbatim
# >          INCY is INTEGER
# >           On entry, INCY specifies the increment for the elements of
# >           CY. INCY must not be zero.
# > \endverbatim
# >
# > \param[in] C
# > \verbatim
# >          C is DOUBLE PRECISION
# >           On entry, C specifies the cosine, cos.
# > \endverbatim
# >
# > \param[in] S
# > \verbatim
# >          S is DOUBLE PRECISION
# >           On entry, S specifies the sine, sin.
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
# > \ingroup complex16_blas_level1
#
#  =====================================================================
from ..util import slice_


def zdrot(N, ZX, INCX, ZY, INCY, C, S):
    #
    #  -- Reference BLAS level1 routine (version 3.7.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     December 2016
    #
    #     .. Scalar Arguments ..
    # INTEGER            INCX, INCY, N
    # DOUBLE PRECISION   C, S
    #     ..
    #     .. Array Arguments ..
    # COMPLEX*16         CX( * ), CY( * )
    #     ..
    #
    # =====================================================================
    if N <= 0:
        return
    x_slice = slice_(N, INCX)
    y_slice = slice_(N, INCY)
    X_TEMP = C * ZX[x_slice] + S * ZY[y_slice]
    ZY[y_slice] = -S * ZX[x_slice] + C * ZY[y_slice]
    ZX[x_slice] = X_TEMP
