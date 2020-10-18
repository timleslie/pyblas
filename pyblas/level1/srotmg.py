# > \brief \b SROTMG
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def SROTMG(SD1,SD2,SX1,SY1,SPARAM)
#
#       .. Scalar Arguments ..
#       REAL SD1,SD2,SX1,SY1
#       ..
#       .. Array Arguments ..
#       REAL SPARAM(5)
#       ..
#
#
# > \par Purpose:
#  =============
# >
# > \verbatim
# >
# >    CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H WHICH ZEROS
# >    THE SECOND COMPONENT OF THE 2-VECTOR  (SQRT(SD1)*SX1,SQRT(SD2)*SY2)**T.
# >    WITH SPARAM(1)=SFLAG, H HAS ONE OF THE FOLLOWING FORMS..
# >
# >    SFLAG=-1.E0     SFLAG=0.E0        SFLAG=1.E0     SFLAG=-2.E0
# >
# >      (SH11  SH12)    (1.E0  SH12)    (SH11  1.E0)    (1.E0  0.E0)
# >    H=(          )    (          )    (          )    (          )
# >      (SH21  SH22),   (SH21  1.E0),   (-1.E0 SH22),   (0.E0  1.E0).
# >    LOCATIONS 2-4 OF SPARAM CONTAIN SH11,SH21,SH12, AND SH22
# >    RESPECTIVELY. (VALUES OF 1.E0, -1.E0, OR 0.E0 IMPLIED BY THE
# >    VALUE OF SPARAM(1) ARE NOT STORED IN SPARAM.)
# >
# >    THE VALUES OF GAMSQ AND RGAMSQ SET IN THE DATA STATEMENT MAY BE
# >    INEXACT.  THIS IS OK AS THEY ARE ONLY USED FOR TESTING THE SIZE
# >    OF SD1 AND SD2.  ALL ACTUAL SCALING OF DATA IS DONE USING GAM.
# >
# > \endverbatim
#
#  Arguments:
#  ==========
#
# > \param[in,out] SD1
# > \verbatim
# >          SD1 is REAL
# > \endverbatim
# >
# > \param[in,out] SD2
# > \verbatim
# >          SD2 is REAL
# > \endverbatim
# >
# > \param[in,out] SX1
# > \verbatim
# >          SX1 is REAL
# > \endverbatim
# >
# > \param[in] SY1
# > \verbatim
# >          SY1 is REAL
# > \endverbatim
# >
# > \param[out] SPARAM
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
def SROTMG(SD1, SD2, SX1, SY1, SPARAM):
    #
    #  -- Reference BLAS level1 routine (version 3.8.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     November 2017
    #
    #     .. Scalar Arguments ..
    # REAL SD1,SD2,SX1,SY1
    #     ..
    #     .. Array Arguments ..
    # REAL SPARAM(5)
    #     ..
    #
    #  =====================================================================
    #
    #     .. Local Scalars ..
    #    REAL GAM,GAMSQ,ONE,RGAMSQ,SFLAG,SH11,SH12,SH21,SH22,SP1,SP2,SQ1,
    #   $     SQ2,STEMP,SU,TWO,ZERO
    #     ..
    #     .. Intrinsic Functions ..
    # INTRINSIC ABS
    #     ..
    #     .. Data statements ..
    #
    # DATA ZERO,ONE,TWO/0.E0,1.E0,2.E0/
    # DATA GAM,GAMSQ,RGAMSQ/4096.E0,1.67772E7,5.96046E-8/
    #     ..

    IDENTITY_MATRIX = -2
    FULL_MATRIX = -1
    IDENT_DIAG = 0
    DIAG_VALUES = 1

    if SD1 < 0:
        SD1 = 0  # H, D?
        SD2 = 0
        SX1 = 0
        SPARAM[:] = [FULL_MATRIX, 0, 0, 0, 0]
        return

    if (
        SD2 * SY1 == 0
    ):  # The vector is already zero, so we can just use the identity matrix!
        SPARAM[:] = [IDENTITY_MATRIX, 0, 0, 0, 0]
        return

    # REGULAR-CASE..
    SQ2 = SD2 * SY1 * SY1
    SQ1 = SD1 * SX1 * SX1
    #
    if abs(SQ1) > abs(SQ2):
        SH21 = -SY1 / SX1
        SH12 = (SD2 * SY1) / (SD1 * SX1)
        SU = 1 - SH12 * SH21
        if SU > 0:
            SFLAG = IDENT_DIAG
            SD1 = SD1 / SU
            SD2 = SD2 / SU
            SX1 = SX1 * SU
        else:
            SPARAM[:] = FULL_MATRIX, 0, 0, 0, 0
            SD1 = 0
            SD2 = 0
            SX1 = 0
            return
    else:
        if SD2 >= 0:
            SFLAG = DIAG_VALUES
            SH11 = (SD1 * SX1) / (SD2 * SY1)
            SH22 = SX1 / SY1
            SU = 1 + SH11 * SH22
            STEMP = SD2 / SU
            SD2 = SD1 / SU
            SD1 = STEMP
            SX1 = SY1 * SU
        else:
            SPARAM[:] = FULL_MATRIX, 0, 0, 0, 0
            SD1 = 0
            SD2 = 0
            SX1 = 0
            return
    #     PROCESURE..SCALE-CHECK
    GAM = 4096
    GAMSQ = 1.67772e7
    RGAMSQ = 5.96046e-8
    if SD1 != 0:
        while (SD1 <= RGAMSQ) or (SD1 >= GAMSQ):
            # First time round we convert to a FULL_MATRIX to allow scaling
            if SFLAG == IDENT_DIAG:
                SH11 = 1
                SH22 = 1
            else:
                SH21 = -1
                SH12 = 1
            SFLAG = FULL_MATRIX

            if SD1 <= RGAMSQ:
                SD1 = SD1 * GAM ** 2
                SX1 = SX1 / GAM
                SH11 = SH11 / GAM
                SH12 = SH12 / GAM
            else:
                SD1 = SD1 / GAM ** 2
                SX1 = SX1 * GAM
                SH11 = SH11 * GAM
                SH12 = SH12 * GAM

    if SD2 != 0:
        while (abs(SD2) <= RGAMSQ) or (abs(SD2) >= GAMSQ):
            # First time round we convert to a FULL_MATRIX to allow scaling
            if SFLAG == IDENT_DIAG:
                SH11 = 1
                SH22 = 1
            else:
                SH21 = -1
                SH12 = 1
            SFLAG = FULL_MATRIX

            if abs(SD2) <= RGAMSQ:
                SD2 = SD2 * GAM ** 2
                SH21 = SH21 / GAM
                SH22 = SH22 / GAM
            else:
                SD2 = SD2 / GAM ** 2
                SH21 = SH21 * GAM
                SH22 = SH22 * GAM

    if SFLAG < 0:
        SPARAM[2] = SH11
        SPARAM[3] = SH21
        SPARAM[4] = SH12
        SPARAM[5] = SH22
    elif SFLAG == IDENT_DIAG:
        SPARAM[3] = SH21
        SPARAM[4] = SH12
    else:
        SPARAM[2] = SH11
        SPARAM[5] = SH22

    SPARAM[1] = SFLAG
