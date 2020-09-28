# > \brief \b DNRM2
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def DNRM2(N,X,INCX)
#
#       .. Scalar Arguments ..
#       INTEGER INCX,N
#       ..
#       .. Array Arguments ..
#       DOUBLE PRECISION X(*)
#       ..
#
#
# > \par Purpose:
#  =============
# >
# > \verbatim
# >
# > DNRM2 returns the euclidean norm of a vector via the function
# > name, so that
# >
# >    DNRM2 := sqrt( x'*x )
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
# > \param[in] X
# > \verbatim
# >          X is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
# > \endverbatim
# >
# > \param[in] INCX
# > \verbatim
# >          INCX is INTEGER
# >         storage spacing between elements of DX
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
# >  -- This version written on 25-October-1982.
# >     Modified on 14-October-1993 to inline the call to DLASSQ.
# >     Sven Hammarling, Nag Ltd.
# > \endverbatim
# >
#  =====================================================================
from math import sqrt


def DNRM2(N, X, INCX):
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
    # DOUBLE PRECISION X(*)
    #     ..
    #
    #  =====================================================================

    if N < 1 or INCX < 1:
        NORM = 0
    elif N == 1:
        NORM = abs(X(1))
    else:
        SCALE = 0
        SSQ = 1
        # The following loop is equivalent to this call to the LAPACK
        # auxiliary routine:
        # CALL DLASSQ( N, X, INCX, SCALE, SSQ )
        for IX in range(0, 1 + (N - 1) * INCX, INCX):
            if X[IX] != 0:
                ABSXI = abs(X[IX])
                if SCALE < ABSXI:
                    SSQ = 1 + SSQ * (SCALE / ABSXI) ** 2
                    SCALE = ABSXI
                else:
                    SSQ = SSQ + (ABSXI / SCALE) ** 2
        NORM = SCALE * sqrt(SSQ)
    return NORM
