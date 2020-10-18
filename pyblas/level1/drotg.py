# > \brief \b DROTG
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def DROTG(DA,DB,C,S)
#
#       .. Scalar Arguments ..
#       DOUBLE PRECISION C,DA,DB,S
#       ..
#
#
# > \par Purpose:
#  =============
# >
# > \verbatim
# >
# >    DROTG construct givens plane rotation.
# > \endverbatim
#
#  Arguments:
#  ==========
#
# > \param[in] DA
# > \verbatim
# >          DA is DOUBLE PRECISION
# > \endverbatim
# >
# > \param[in] DB
# > \verbatim
# >          DB is DOUBLE PRECISION
# > \endverbatim
# >
# > \param[out] C
# > \verbatim
# >          C is DOUBLE PRECISION
# > \endverbatim
# >
# > \param[out] S
# > \verbatim
# >          S is DOUBLE PRECISION
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


def drotg(DA, DB):
    #
    #  -- Reference BLAS level1 routine (version 3.8.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     November 2017
    #
    #     .. Scalar Arguments ..
    # DOUBLE PRECISION C,DA,DB,S
    #     ..
    #
    #  =====================================================================
    #
    #     .. Local Scalars ..
    # DOUBLE PRECISION R,ROE,SCALE,Z
    #     ..
    #     .. Intrinsic Functions ..
    # INTRINSIC DABS,DSIGN,DSQRT
    #     ..
    SCALE = abs(DA) + abs(DB)
    if SCALE == 0:
        C, S, R, Z = 1.0, 0.0, 0.0, 0.0
    else:
        ROE = DA if abs(DA) > abs(DB) else DB
        R = np.sign(ROE) * SCALE * sqrt((DA / SCALE) ** 2 + (DB / SCALE) ** 2)
        C = DA / R
        S = DB / R
        if abs(DA) > abs(DB):
            Z = S
        elif abs(DB) >= abs(DA) and C != 0.0:
            Z = 1.0 / C
        else:
            Z = 1.0
    return C, S, R, Z
