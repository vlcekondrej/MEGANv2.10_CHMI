!-----------------------------------------------------------------------
!   SUBROUTINE: SOLARANGLE
!
!   Description: To calculate the solar zenith angle.  This will give
!                SIN(BETA), not the BETA.
!
!   Call: None
!
!   Require: None
!
!   Input:
!            1) Day of year
!            2) Latitude
!            3) Hour
!
!   Output: CALCBETA (Solar zenith angle)
!
!   Created by Tan 11/15/06  (based on xxxx's program)
!
!-----------------------------------------------------------------------
      SUBROUTINE SOLARANGLE( DAY, SHOUR, LAT, SINbeta)

      IMPLICIT NONE

! input
      INTEGER,INTENT(IN) :: DAY     ! DOY or julian day
      REAL,INTENT(IN) :: SHOUR      ! Solar hour
      REAL,INTENT(IN) :: LAT        ! Latitude
! output
      REAL,INTENT(OUT) :: SINbeta
! Local
      REAL    :: BETA                 ! solar elevation angle
      REAL    :: sindelta, cosdelta, A, B
! Constants
      REAL,PARAMETER :: PI    = 3.14159,  
     &                  D2RAD = PI/180.0, 
     &                  RAD2D = 180.0/PI

! Calculation
      sindelta = -SIN(0.40907) * COS( 6.28*(DAY+10)/365 )
      cosdelta = (1-sindelta**2.)**0.5

      A = SIN( LAT*D2RAD ) * sindelta
      B = COS( LAT*D2RAD ) * cosdelta

      SINbeta = A + B * COS(2*PI*(SHOUR-12)/24)  ! This will be transfered
                                                 ! to gamma_p function

      BETA = ASIN(sinbeta)*RAD2D    ! This is not used.

      RETURN
      END SUBROUTINE SOLARANGLE
!-----------------------------------------------------------------------

