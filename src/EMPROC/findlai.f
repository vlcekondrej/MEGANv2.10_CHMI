!-----------------------------------------------------------------------
!   SUBROUTINE: FINDLAI
!
!   Description: Find current LAI and previous LAI from LAIS46 and IDATE
!
!   Call: None
!
!   Require: None
!
!   Input:
!            1) IDATE : current julian date
!
!   Output:  1) LAIp_I
!            2) LAIc_I
!
!   Created by Tan 07/28/11
!
!-----------------------------------------------------------------------
      SUBROUTINE FINDLAI( IDATE, MXLAI, LAIp_DY, LAIp_HR, LAIc_DY, LAIc_HR)

      IMPLICIT NONE

! input
      INTEGER,INTENT(IN) :: IDATE  ! YYYYJJJ
      INTEGER,INTENT(IN) :: MXLAI
! output
      INTEGER,INTENT(OUT) :: LAIp_DY, LAIp_HR
      INTEGER,INTENT(OUT) :: LAIc_DY, LAIc_HR
! Local
      INTEGER :: JJJ
      INTEGER :: LAIp_I, LAIc_I
      REAL    :: XXX

! Calculation

      JJJ = MOD(IDATE,1000)
      XXX = JJJ/8.0
      LAIc_I = CEILING(XXX)

      IF (LAIc_I .EQ. 1) THEN
        LAIp_I = MXLAI
      ELSE
        LAIp_I = LAIc_I - 1
      ENDIF

      LAIc_HR = (MOD(LAIc_I-1,24))*10000
      XXX     = (LAIc_I-1)/24.0
      LAIc_DY = FLOOR(XXX)
      LAIp_HR = (MOD(LAIp_I-1,24))*10000
      XXX     = (LAIp_I-1)/24.0
      LAIp_DY = FLOOR(XXX)

      RETURN
      END SUBROUTINE FINDLAI
!-----------------------------------------------------------------------
