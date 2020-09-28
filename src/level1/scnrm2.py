# > \brief \b SCNRM2
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def SCNRM2(N,X,INCX)
#
#       .. Scalar Arguments ..
#       INTEGER INCX,N
#       ..
#       .. Array Arguments ..
#       COMPLEX X(*)
#       ..
#
#
# > \par Purpose:
#  =============
# >
# > \verbatim
# >
# > SCNRM2 returns the euclidean norm of a vector via the function
# > name, so that
# >
# >    SCNRM2 := sqrt( x**H*x )
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
# >          X is COMPLEX array, dimension (N)
# >         complex vector with N elements
# > \endverbatim
# >
# > \param[in] INCX
# > \verbatim
# >          INCX is INTEGER
# >         storage spacing between elements of X
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
# >  -- This version written on 25-October-1982.
# >     Modified on 14-October-1993 to inline the call to CLASSQ.
# >     Sven Hammarling, Nag Ltd.
# > \endverbatim
# >
#  =====================================================================
from math import sqrt


def SCNRM2(N, X, INCX):
    #
    #  -- Reference BLAS level1 routine (version 3.8.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     November 2017
    #
    #     .. Scalar Arguments ..
    #   INTEGER INCX,N
    #     ..
    #     .. Array Arguments ..
    #   COMPLEX X(*)
    #     ..
    #
    #  =====================================================================
    #
    #     .. Parameters ..
    #   REAL ONE,ZERO
    #   PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
    #     ..
    #     .. Local Scalars ..
    #   REAL NORM,SCALE,SSQ,TEMP
    #   INTEGER IX
    #     ..
    #     .. Intrinsic Functions ..
    #   INTRINSIC ABS,AIMAG,REAL,SQRT
    #     ..
    if N < 1 or INCX < 1:
        NORM = 0
    else:
        SCALE = 0
        SSQ = 1
        #        The following loop is equivalent to this call to the LAPACK
        #        auxiliary routine:
        #        CALL CLASSQ( N, X, INCX, SCALE, SSQ )
        #
        for IX in range(0, 1 + (N - 1) * INCX, INCX):
            if (X[IX]).real != 0:
                TEMP = abs((X[IX]).real)
                if SCALE < TEMP:
                    SSQ = 1 + SSQ * (SCALE / TEMP) ** 2
                    SCALE = TEMP
                else:
                    SSQ = SSQ + (TEMP / SCALE) ** 2
            if (X[IX]).image != 0:
                TEMP = abs((X[IX]).imag)
                if SCALE < TEMP:
                    SSQ = 1 + SSQ * (SCALE / TEMP) ** 2
                    SCALE = TEMP
                else:
                    SSQ = SSQ + (TEMP / SCALE) ** 2
        NORM = SCALE * sqrt(SSQ)
    #
    return NORM
