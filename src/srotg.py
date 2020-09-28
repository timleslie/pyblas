# > \brief \b SROTG
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def SROTG(SA,SB,C,S)
#
#       .. Scalar Arguments ..
#       REAL C,S,SA,SB
#       ..
#
#
# > \par Purpose:
#  =============
# >
# > \verbatim
# >
# >    SROTG construct givens plane rotation.
# > \endverbatim
#
#  Arguments:
#  ==========
#
# > \param[in] SA
# > \verbatim
# >          SA is REAL
# > \endverbatim
# >
# > \param[in] SB
# > \verbatim
# >          SB is REAL
# > \endverbatim
# >
# > \param[out] C
# > \verbatim
# >          C is REAL
# > \endverbatim
# >
# > \param[out] S
# > \verbatim
# >          S is REAL
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
# > \endverbatim
# >
#  =====================================================================
from math import sqrt
from util import sign


def SROTG(SA, SB, C, S):
    #
    #  -- Reference BLAS level1 routine (version 3.8.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     November 2017
    #
    #     .. Scalar Arguments ..
    # REAL C,S,SA,SB
    #     ..
    #
    #  =====================================================================
    #
    #     .. Local Scalars ..
    # REAL R,ROE,SCALE,Z
    #     ..
    #     .. Intrinsic Functions ..
    # INTRINSIC ABS,SIGN,SQRT
    #     ..
    ROE = SB
    if abs(SA) > abs(SB):
        ROE = SA
    SCALE = abs(SA) + abs(SB)
    if SCALE == 0.0:
        C = 1.0
        S = 0.0
        R = 0.0
        Z = 0.0
    else:
        R = SCALE * sqrt((SA / SCALE) ** 2 + (SB / SCALE) ** 2)
        R = sign(ROE) * R
        C = SA / R
        S = SB / R
        Z = 1.0
        if abs(SA) > abs(SB):
            Z = S
        if abs(SB) >= abs(SA) and C != 0.0:
            Z = 1.0 / C
    SA = R
    SB = Z
