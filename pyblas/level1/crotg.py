# > \brief \b CROTG
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def CROTG(CA,CB,C,S)
#
#       .. Scalar Arguments ..
#       COMPLEX CA,CB,S
#       REAL C
#       ..
#
#
# > \par Purpose:
#  =============
# >
# > \verbatim
# >
# > CROTG determines a complex Givens rotation.
# > \endverbatim
#
#  Arguments:
#  ==========
#
# > \param[in] CA
# > \verbatim
# >          CA is COMPLEX
# > \endverbatim
# >
# > \param[in] CB
# > \verbatim
# >          CB is COMPLEX
# > \endverbatim
# >
# > \param[out] C
# > \verbatim
# >          C is REAL
# > \endverbatim
# >
# > \param[out] S
# > \verbatim
# >          S is COMPLEX
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
# > \ingroup complex_blas_level1
#
#  =====================================================================
from math import sqrt


def CROTG(CA, CB, C, S):
    #
    #  -- Reference BLAS level1 routine (version 3.8.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     November 2017
    #
    #     .. Scalar Arguments ..
    # COMPLEX CA,CB,S
    # REAL C
    #     ..
    #
    #  =====================================================================
    #
    #     .. Local Scalars ..
    # COMPLEX ALPHA
    # REAL NORM,SCALE
    #     ..
    #     .. Intrinsic Functions ..
    # INTRINSIC CABS,CONJG,SQRT
    #     ..
    if abs(CA) == 0.0:
        C = 0.0
        S = 1 + 0j
        CA = CB
    else:
        SCALE = abs(CA) + abs(CB)
        NORM = SCALE * sqrt((abs(CA / SCALE)) ** 2 + (abs(CB / SCALE)) ** 2)
        ALPHA = CA / abs(CA)
        C = abs(CA) / NORM
        S = ALPHA * (CB).conjugate() / NORM
        CA = ALPHA * NORM
