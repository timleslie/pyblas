# > \brief \b DROTM
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def DROTM(N,DX,INCX,DY,INCY,DPARAM)
#
#       .. Scalar Arguments ..
#       INTEGER INCX,INCY,N
#       ..
#       .. Array Arguments ..
#       DOUBLE PRECISION DPARAM(5),DX(*),DY(*)
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
# >    (DX**T) , WHERE **T INDICATES TRANSPOSE. THE ELEMENTS OF DX ARE IN
# >    (DY**T)
# >
# >    DX(LX+I*INCX), I = 0 TO N-1, WHERE LX = 1 IF INCX >= 0, else:
# >    LX = (-INCX)*N, AND SIMILARLY FOR SY USING LY AND INCY.
# >    WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS..
# >
# >    DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0
# >
# >      (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0)
# >    H=(          )    (          )    (          )    (          )
# >      (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0).
# >    SEE DROTMG FOR A DESCRIPTION OF DATA STORAGE IN DPARAM.
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
# > \param[in,out] DX
# > \verbatim
# >          DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
# > \endverbatim
# >
# > \param[in] INCX
# > \verbatim
# >          INCX is INTEGER
# >         storage spacing between elements of DX
# > \endverbatim
# >
# > \param[in,out] DY
# > \verbatim
# >          DY is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
# > \endverbatim
# >
# > \param[in] INCY
# > \verbatim
# >          INCY is INTEGER
# >         storage spacing between elements of DY
# > \endverbatim
# >
# > \param[in] DPARAM
# > \verbatim
# >          DPARAM is DOUBLE PRECISION array, dimension (5)
# >     DPARAM(1)=DFLAG
# >     DPARAM(2)=DH11
# >     DPARAM(3)=DH21
# >     DPARAM(4)=DH12
# >     DPARAM(5)=DH22
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
# > \ingroup double_blas_level1
#
#  =====================================================================
from ..util import slice_


def drotm(N, DX, INCX, DY, INCY, DPARAM):
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
    # DOUBLE PRECISION DPARAM(5),DX(*),DY(*)
    #     ..
    #
    #  =====================================================================
    #
    #     .. Local Scalars ..
    # DOUBLE PRECISION DFLAG,DH11,DH12,DH21,DH22,TWO,W,Z,ZERO
    # INTEGER I,KX,KY,NSTEPS
    #     ..
    #     .. Data statements ..
    # DATA ZERO,TWO/0.D0,2.D0/
    #     ..
    #
    [DFLAG, DH11, DH21, DH12, DH22] = DPARAM
    if N <= 0:
        return
    if DFLAG == -2:  # Identity matrix
        return

    sx = slice_(N, INCX)
    sy = slice_(N, INCY)
    if DFLAG == -1:  # Full matrix
        TMP_X = DH11 * DX[sx] + DH12 * DY[sy]
        DY[sy] = DH21 * DX[sx] + DH22 * DY[sy]
        DX[sx] = TMP_X
    elif DFLAG == 0:  # Off-diagonal
        TMP_X = DX[sx] + DH12 * DY[sy]
        DY[sy] = DH21 * DX[sx] + DY[sy]
        DX[sx] = TMP_X
    elif DFLAG == 1:
        TMP_X = DH11 * DX[sx] + DY[sy]
        DY[sy] = -DX[sx] + DH22 * DY[sy]
        DX[sx] = TMP_X
