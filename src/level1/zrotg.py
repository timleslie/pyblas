# > \brief \b ZROTG
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def ZROTG(CA,CB,C,S)
#
#       .. Scalar Arguments ..
#       COMPLEX*16 CA,CB,S
#       DOUBLE PRECISION C
#       ..
#
#
# > \par Purpose:
#  =============
# >
# > \verbatim
# >
# >    ZROTG determines a double complex Givens rotation.
# > \endverbatim
#
#  Arguments:
#  ==========
#
# > \param[in] CA
# > \verbatim
# >          CA is COMPLEX*16
# > \endverbatim
# >
# > \param[in] CB
# > \verbatim
# >          CB is COMPLEX*16
# > \endverbatim
# >
# > \param[out] C
# > \verbatim
# >          C is DOUBLE PRECISION
# > \endverbatim
# >
# > \param[out] S
# > \verbatim
# >          S is COMPLEX*16
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
# > \ingroup complex16_blas_level1
#
#  =====================================================================
from math import sqrt


def ZROTG(CA, CB, C, S):
    #
    #  -- Reference BLAS level1 routine (version 3.8.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     November 2017
    #
    #     .. Scalar Arguments ..
    # COMPLEX*16 CA,CB,S
    # DOUBLE PRECISION C
    #     ..
    #
    #  =====================================================================
    #
    #     .. Local Scalars ..
    # COMPLEX*16 ALPHA
    # DOUBLE PRECISION NORM,SCALE
    #     ..
    #     .. Intrinsic Functions ..
    # INTRINSIC CDABS,DCMPLX,DCONJG,DSQRT
    #     ..
    if abs(CA) == 0:
        C = 0
        S = 0
        CA = CB
    else:
        SCALE = abs(CA) + abs(CB)
        NORM = SCALE * sqrt((abs(CA / SCALE)) ** 2 + (abs(CB / SCALE)) ** 2)
        ALPHA = CA / abs(CA)
        C = abs(CA) / NORM
        S = ALPHA * (CB).conjugate() / NORM
        CA = ALPHA * NORM
