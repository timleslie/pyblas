# > \brief \b SROTM
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def SROTM(N,SX,INCX,SY,INCY,SPARAM)
#
#       .. Scalar Arguments ..
#       INTEGER INCX,INCY,N
#       ..
#       .. Array Arguments ..
#       REAL SPARAM(5),SX(*),SY(*)
#       ..
#
#
# > \par Purpose:
#  =============
# >
# > \verbatim
# >
# >    APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N MATRIX
# >
# >    (SX**T) , WHERE **T INDICATES TRANSPOSE. THE ELEMENTS OF SX ARE IN
# >    (SX**T)
# >
# >    SX(LX+I*INCX), I = 0 TO N-1, WHERE LX = 1 IF INCX >= 0, else:
# >    LX = (-INCX)*N, AND SIMILARLY FOR SY USING USING LY AND INCY.
# >    WITH SPARAM(1)=SFLAG, H HAS ONE OF THE FOLLOWING FORMS..
# >
# >    SFLAG=-1.E0     SFLAG=0.E0        SFLAG=1.E0     SFLAG=-2.E0
# >
# >      (SH11  SH12)    (1.E0  SH12)    (SH11  1.E0)    (1.E0  0.E0)
# >    H=(          )    (          )    (          )    (          )
# >      (SH21  SH22),   (SH21  1.E0),   (-1.E0 SH22),   (0.E0  1.E0).
# >    SEE  SROTMG FOR A DESCRIPTION OF DATA STORAGE IN SPARAM.
# >
# > \endverbatim
#
#  Arguments:
#  ==========
#
# > \param[in] N
# > \verbatim
# >          N is INTEGER
# >         number of elements in input vector(s)
# > \endverbatim
# >
# > \param[in,out] SX
# > \verbatim
# >          SX is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
# > \endverbatim
# >
# > \param[in] INCX
# > \verbatim
# >          INCX is INTEGER
# >         storage spacing between elements of SX
# > \endverbatim
# >
# > \param[in,out] SY
# > \verbatim
# >          SY is REAL array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
# > \endverbatim
# >
# > \param[in] INCY
# > \verbatim
# >          INCY is INTEGER
# >         storage spacing between elements of SY
# > \endverbatim
# >
# > \param[in] SPARAM
# > \verbatim
# >          SPARAM is REAL array, dimension (5)
# >     SPARAM(1)=SFLAG
# >     SPARAM(2)=SH11
# >     SPARAM(3)=SH21
# >     SPARAM(4)=SH12
# >     SPARAM(5)=SH22
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
# > \date November 2017
#
# > \ingroup single_blas_level1
#
#  =====================================================================
from ..util import slice_


def srotm(N, SX, INCX, SY, INCY, SPARAM):
    #
    #  -- Reference BLAS level1 routine (version 3.8.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     November 2017
    #
    #     .. Scalar Arguments ..
    # INTEGER INCX,INCY,N
    #     ..
    #     .. Array Arguments ..
    # REAL SPARAM(5),SX(*),SY(*)
    #     ..
    #
    #  =====================================================================
    #
    #     .. Local Scalars ..
    # REAL SFLAG,SH11,SH12,SH21,SH22,TWO,W,Z,ZERO
    # INTEGER I,KX,KY,NSTEPS
    #     ..
    #     .. Data statements ..
    # DATA ZERO,TWO/0.E0,2.E0/
    #     ..
    #
    [SFLAG, SH11, SH21, SH12, SH22] = SPARAM
    if N <= 0:
        return
    if SFLAG == -2:  # Identity matrix
        return

    sx = slice_(N, INCX)
    sy = slice_(N, INCY)
    if SFLAG == -1:  # Full matrix
        TMP_X = SH11 * SX[sx] + SH12 * SY[sy]
        SY[sy] = SH21 * SX[sx] + SH22 * SY[sy]
        SX[sx] = TMP_X
    elif SFLAG == 0:  # Off-diagonal
        TMP_X = SX[sx] + SH12 * SY[sy]
        SY[sy] = SH21 * SX[sx] + SY[sy]
        SX[sx] = TMP_X
    elif SFLAG == 1:
        TMP_X = SH11 * SX[sx] + SY[sy]
        SY[sy] = -SX[sx] + SH22 * SY[sy]
        SX[sx] = TMP_X
