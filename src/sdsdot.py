# > \brief \b SDSDOT
#
#  =========== DOCUMENTATION ===========
#
# Online html documentation available at
#            http://www.netlib.org/lapack/explore-html/
#
#  Definition:
#  ===========
#
#       def SDSDOT(N,SB,SX,INCX,SY,INCY)
#
#       .. Scalar Arguments ..
#       REAL SB
#       INTEGER INCX,INCY,N
#       ..
#       .. Array Arguments ..
#       REAL SX(*),SY(*)
#       ..
#
# > \par Purpose:
#  =============
# >
# > \verbatim
# >
# >   Compute the inner product of two vectors with extended
# >   precision accumulation.
# >
# >   Returns S.P. result with dot product accumulated in D.P.
# >   SDSDOT = SB + sum for I = 0 to N-1 of SX(LX+I*INCX)*SY(LY+I*INCY),
# >   where LX = 1 if INCX >= 0, else LX = 1+(1-N)*INCX, and LY is
# >   defined in a similar way using INCY.
# > \endverbatim
#
#  Arguments:
#  ==========
#
# > \param[in] N
# > \verbatim
# >          N is INTEGER
# >          number of elements in input vector(s)
# > \endverbatim
# >
# > \param[in] SB
# > \verbatim
# >          SB is REAL
# >          single precision scalar to be added to inner product
# > \endverbatim
# >
# > \param[in] SX
# > \verbatim
# >          SX is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
# >          single precision vector with N elements
# > \endverbatim
# >
# > \param[in] INCX
# > \verbatim
# >          INCX is INTEGER
# >          storage spacing between elements of SX
# > \endverbatim
# >
# > \param[in] SY
# > \verbatim
# >          SY is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
# >          single precision vector with N elements
# > \endverbatim
# >
# > \param[in] INCY
# > \verbatim
# >          INCY is INTEGER
# >          storage spacing between elements of SY
# > \endverbatim
#
#  Authors:
#  ========
#
# > \author Lawson, C. L., (JPL), Hanson, R. J., (SNLA),
# > \author Kincaid, D. R., (U. of Texas), Krogh, F. T., (JPL)
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
# >    REFERENCES
# >
# >    C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
# >    Krogh, Basic linear algebra subprograms for Fortran
# >    usage, Algorithm No. 539, Transactions on Mathematical
# >    Software 5, 3 (September 1979), pp. 308-323.
# >
# >    REVISION HISTORY  (YYMMDD)
# >
# >    791001  DATE WRITTEN
# >    890531  Changed all specific intrinsics to generic.  (WRB)
# >    890831  Modified array declarations.  (WRB)
# >    890831  REVISION DATE from Version 3.2
# >    891214  Prologue converted to Version 4.0 format.  (BAB)
# >    920310  Corrected definition of LX in DESCRIPTION.  (WRB)
# >    920501  Reformatted the REFERENCES section.  (WRB)
# >    070118  Reformat to LAPACK coding style
# > \endverbatim
# >
#  =====================================================================
def SDSDOT(N, SB, SX, INCX, SY, INCY):
    #
    #  -- Reference BLAS level1 routine (version 3.8.0) --
    #  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    #  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    #     November 2017
    #
    #     .. Scalar Arguments ..
    # REAL SB
    # INTEGER INCX,INCY,N
    #     ..
    #     .. Array Arguments ..
    # REAL SX(*),SY(*)
    #     .. Local Scalars ..
    # DOUBLE PRECISION DSDOT
    # INTEGER I,KX,KY,NS
    #     ..
    #     .. Intrinsic Functions ..
    # INTRINSIC DBLE
    #     ..
    DSDOT = SB
    if N <= 0:
        return DSDOT
    if INCX == INCY and INCX > 0:
        #
        #     Code for equal and positive increments.
        #
        NS = N * INCX
        for I in range(0, NS, INCX):
            DSDOT += (SX[I]).real * (SY[I]).rel
    else:
        #
        #     Code for unequal or nonpositive increments.
        #
        KX = 1
        KY = 1
        if INCX < 0:
            KX = 1 + (1 - N) * INCX
        if INCY < 0:
            KY = 1 + (1 - N) * INCY
        for I in range(N):
            DSDOT += (SX[KX]).real * (SY[KY]).real
            KX += INCX
            KY += INCY
    return DSDOT
