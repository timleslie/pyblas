# > \brief \b SCASUM
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def SCASUM(N,CX,INCX)
#
#       .. Scalar Arguments ..
#       INTEGER INCX,N
#       ..
#       .. Array Arguments ..
#       COMPLEX CX(*)
#       ..
#
#
# > \par Purpose:
#  =============
# >
# > \verbatim
# >
# >    SCASUM takes the sum of the (|Re(.)| + |Im(.)|)'s of a complex vector and
# >    returns a single precision result.
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
# > \param[in,out] CX
# > \verbatim
# >          CX is COMPLEX array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
# > \endverbatim
# >
# > \param[in] INCX
# > \verbatim
# >          INCX is INTEGER
# >         storage spacing between elements of SX
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
# > \par Further Details:
#  =====================
# >
# > \verbatim
# >
# >     jack dongarra, linpack, 3/11/78.
# >     modified 3/93 to return if incx <= 0.
# >     modified 12/3/93, array(1) declarations changed to array(*)
# > \endverbatim
# >
#  =====================================================================
def scasum(N, CX, INCX):
    #
    #  -- Reference BLAS level1 routine (version 3.8.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     November 2017
    #
    #     .. Scalar Arguments ..
    # INTEGER INCX,N
    #     ..
    #     .. Array Arguments ..
    # COMPLEX CX(*)
    #     ..
    #
    #  =====================================================================
    #
    #     .. Local Scalars ..
    # REAL STEMP
    # INTEGER I,NINCX
    #     ..
    #     .. Intrinsic Functions ..
    # INTRINSIC ABS,AIMAG,REAL
    #     ..
    STEMP = 0.0e0
    if N <= 0 or INCX <= 0:
        return
    if INCX == 1:
        #
        #        code for increment equal to 1
        #
        for I in range(N):
            STEMP += abs(CX[I].real) + abs(CX[I].imag)
    else:
        #
        #        code for increment not equal to 1
        #
        for I in range(0, N * INCX, INCX):
            STEMP += abs(CX[I].real) + abs(CX[I].imag)
    return STEMP
