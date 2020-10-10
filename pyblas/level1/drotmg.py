# > \brief \b DROTMG
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def DROTMG(DD1,DD2,DX1,DY1,DPARAM)
#
#       .. Scalar Arguments ..
#       DOUBLE PRECISION DD1,DD2,DX1,DY1
#       ..
#       .. Array Arguments ..
#       DOUBLE PRECISION DPARAM(5)
#       ..
#
#
# > \par Purpose:
#  =============
# >
# > \verbatim
# >
# >    CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H WHICH ZEROS
# >    THE SECOND COMPONENT OF THE 2-VECTOR  (sqrt(DD1)*DX1,sqrt(DD2)*>    DY2)**T.
# >    WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS..
# >
# >    DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0
# >
# >      (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0)
# >    H=(          )    (          )    (          )    (          )
# >      (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0).
# >    LOCATIONS 2-4 OF DPARAM CONTAIN DH11, DH21, DH12, AND DH22
# >    RESPECTIVELY. (VALUES OF 1.D0, -1.D0, OR 0.D0 IMPLIED BY THE
# >    VALUE OF DPARAM(1) ARE NOT STORED IN DPARAM.)
# >
# >    THE VALUES OF GAMSQ AND RGAMSQ SET IN THE DATA STATEMENT MAY BE
# >    INEXACT.  THIS IS OK AS THEY ARE ONLY USED FOR TESTING THE SIZE
# >    OF DD1 AND DD2.  ALL ACTUAL SCALING OF DATA IS DONE USING GAM.
# >
# > \endverbatim
#
#  Arguments:
#  ==========
#
# > \param[in,out] DD1
# > \verbatim
# >          DD1 is DOUBLE PRECISION
# > \endverbatim
# >
# > \param[in,out] DD2
# > \verbatim
# >          DD2 is DOUBLE PRECISION
# > \endverbatim
# >
# > \param[in,out] DX1
# > \verbatim
# >          DX1 is DOUBLE PRECISION
# > \endverbatim
# >
# > \param[in] DY1
# > \verbatim
# >          DY1 is DOUBLE PRECISION
# > \endverbatim
# >
# > \param[out] DPARAM
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
def DROTMG(DD1, DD2, DX1, DY1, DPARAM):
    #
    #  -- Reference BLAS level1 routine (version 3.8.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     November 2017
    #
    #     .. Scalar Arguments ..
    # DOUBLE PRECISION DD1,DD2,DX1,DY1
    #     ..
    #     .. Array Arguments ..
    # DOUBLE PRECISION DPARAM(5)
    #     ..
    #
    #  =====================================================================
    #
    #     .. Local Scalars ..
    #    DOUBLE PRECISION DFLAG,DH11,DH12,DH21,DH22,DP1,DP2,DQ1,DQ2,DTEMP,
    #   $                 DU,GAM,GAMSQ,ONE,RGAMSQ,TWO,ZERO
    #     ..
    #     .. Intrinsic Functions ..
    # INTRINSIC DABS
    #     ..
    #     .. Data statements ..
    #
    # DATA ZERO,ONE,TWO/0.D0,1.D0,2.D0/
    #  DATA GAM,GAMSQ,RGAMSQ/4096.D0,16777216.D0,5.9604645D-8/
    #     ..
    GAM = 4096
    GAMSQ = 16777216
    RGAMSQ = 5.9604645e-8
    if DD1 < 0:
        #        GO ZERO-H-D-AND-DX1..
        DFLAG = -1
        DH11 = 0
        DH12 = 0
        DH21 = 0
        DH22 = 0
        #
        DD1 = 0
        DD2 = 0
        DX1 = 0
    else:
        #        CASE-DD1-NONNEGATIVE
        DP2 = DD2 * DY1
        if DP2 == 0:
            DFLAG = -2
            DPARAM[1] = DFLAG
            return
        #        REGULAR-CASE..
        DP1 = DD1 * DX1
        DQ2 = DP2 * DY1
        DQ1 = DP1 * DX1
        #
        if abs(DQ1) > abs(DQ2):
            DH21 = -DY1 / DX1
            DH12 = DP2 / DP1
            #
            DU = 1 - DH12 * DH21
            #
            if DU > 0:
                DFLAG = 0
                DD1 = DD1 / DU
                DD2 = DD2 / DU
                DX1 = DX1 * DU
        else:

            if DQ2 < 0:
                #               GO ZERO-H-D-AND-DX1..
                DFLAG = -1
                DH11 = 0
                DH12 = 0
                DH21 = 0
                DH22 = 0
                #
                DD1 = 0
                DD2 = 0
                DX1 = 0
            else:
                DFLAG = 1
                DH11 = DP1 / DP2
                DH22 = DX1 / DY1
                DU = 1 + DH11 * DH22
                DTEMP = DD2 / DU
                DD2 = DD1 / DU
                DD1 = DTEMP
                DX1 = DY1 * DU

        # PROCEDURE..SCALE-CHECK
        if DD1 != 0:
            while (DD1 <= RGAMSQ) or (DD1 >= GAMSQ):
                if DFLAG == 0:
                    DH11 = 1
                    DH22 = 1
                    DFLAG = -1
                else:
                    DH21 = -1
                    DH12 = 1
                    DFLAG = -1
                if DD1 <= RGAMSQ:
                    DD1 = DD1 * GAM ** 2
                    DX1 = DX1 / GAM
                    DH11 = DH11 / GAM
                    DH12 = DH12 / GAM
                else:
                    DD1 = DD1 / GAM ** 2
                    DX1 = DX1 * GAM
                    DH11 = DH11 * GAM
                    DH12 = DH12 * GAM

        if DD2 != 0:
            while (abs(DD2) <= RGAMSQ) or (abs(DD2) >= GAMSQ):
                if DFLAG == 0:
                    DH11 = 1
                    DH22 = 1
                    DFLAG = -1
                else:
                    DH21 = -1
                    DH12 = 1
                    DFLAG = -1
                if abs(DD2) <= RGAMSQ:
                    DD2 = DD2 * GAM ** 2
                    DH21 = DH21 / GAM
                    DH22 = DH22 / GAM
                else:
                    DD2 = DD2 / GAM ** 2
                    DH21 = DH21 * GAM
                    DH22 = DH22 * GAM

    if DFLAG < 0:
        DPARAM[2] = DH11
        DPARAM[3] = DH21
        DPARAM[4] = DH12
        DPARAM[5] = DH22
    elif DFLAG == 0:
        DPARAM[3] = DH21
        DPARAM[4] = DH12
    else:
        DPARAM[2] = DH11
        DPARAM[5] = DH22

    DPARAM[1] = DFLAG
