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
import numpy as np


def srotg(SA, SB):
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
    SCALE = abs(SA) + abs(SB)
    if SCALE == 0.0:
        C, S, R, Z = 1.0, 0.0, 0.0, 0.0
    else:
        ROE = SA if abs(SA) > abs(SB) else SB
        R = np.sign(ROE) * SCALE * sqrt((SA / SCALE) ** 2 + (SB / SCALE) ** 2)
        C = SA / R
        S = SB / R
        if abs(SA) > abs(SB):
            Z = S
        elif abs(SB) >= abs(SA) and C != 0.0:
            Z = 1.0 / C
        else:
            Z = 1.0
    return C, S, R, Z
