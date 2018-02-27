!=======================================================================
!     MODULE GAMMA
!
!     This module contain functions to calculate
!     GAMMA_P, GAMMA_T, GAMMA_L, GAMMA_A for BVOCs.
!
!     CONTAINS: 1)GAMMA_LAI
!               2)GAMMA_P
!               3)GAMMA_TLD
!               4)GAMMA_TLI
!               5)GAMMA_A
!               6)GAMMA_S
!               7)GAMMA_CO2
!               8)GAMMA_LAIbidir
!
!     Note:
!
!     Requirement:
!
!     CALLS: SOLARANGLE
!
!     Created by Tan 11/21/06 for MEGAN v2.0
!
!     History:
!     08/01/07 Guenther A. - Move to MEGANv2.02 with modification to
!                            correct calculation of GAMMA_P
!
!=======================================================================

      MODULE GAMMA_etc

      IMPLICIT NONE

!...  Program I/O parameters

!...  External parameters

      CONTAINS
!***********************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Scientific algorithm
!
!             Emission = [EF][GAMMA][RHO]
!           where [EF]    = emission factor (ug/m2h)
!                 [GAMMA] = emission activity factor (non-dimension)
!                 [RHO]   = production and loss within plant canopies
!                           (non-dimensino)
!                 Assumption: [RHO] = 1 (11/27/06) (See PDT_LOT_CP.EXT)
!
!             GAMMA  = [GAMMA_CE][GAMMA_age][GAMMA_SM]
!           where [GAMMA_CE]  = canopy correction factor
!                 [GAMMA_age] = leaf age correction factor
!                 [GAMMA_SM]  = soil moisture correction factor
!                 Assumption: [GAMMA_SM]  = 1 (11/27/06)
!
!             GAMMA_CE = [GAMMA_LAI][GAMMA_P][GAMMA_T]
!           where [GAMMA_LAI] = leaf area index factor
!                 [GAMMA_P]   = PPFD emission activity factor
!                 [GAMMA_T]   = temperature response factor
!
!             Emission = [EF][GAMMA_LAI][GAMMA_P][GAMMA_T][GAMMA_age][GAMMA_SM]
!        Derivation:
!             Emission = [EF][GAMMA_etc](1-LDF) + [EF][GAMMA_etc][LDF][GAMMA_P]
!             Emission = [EF][GAMMA_etc]{ (1-LDF) + [LDF][GAMMA_P] }
!             Emission = [EF][GAMMA_ect]{ (1-LDF) + [LDF][GAMMA_P] }
!           where LDF = light dependent function (non-dimension)
!
!     For ISOPRENE
!                 Assumption: LDF = 1 for isoprene            (11/27/06)
!
!        Final Equation
!             Emission = [EF][GAMMA_LAI][GAMMA_P][GAMMA_T][GAMMA_age][GAMMA_SM]
!
!     For NON-ISOPRENE
!        Final Equation
!             Emission = [EF][GAMMA_LAI][GAMMA_T][GAMMA_age][GAMMA_SM]*
!                        { (1-LDF) + [LDF][GAMMA_P] }
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!=======================================================================
!...  Begin module
!=======================================================================


!-----------------------------------------------------------------------
!.....1) Calculate GAM_L (GAMMA_LAI)
!-----------------------------------------------------------------------
!                            0.49[LAI]
!             GAMMA_LAI = ----------------    (non-dimension)
!                         (1+0.2LAI^2)^0.5
!
!     SUBROUTINE GAMMA_LAI returns the GAMMA_LAI values
!-----------------------------------------------------------------------
      SUBROUTINE GAMMA_LAI( ncols, nrows, lai, gam_l )

      IMPLICIT NONE
! input
      INTEGER,INTENT(IN) :: ncols, nrows
      REAL,DIMENSION(ncols,nrows),INTENT(IN) :: lai
! output
      REAL,DIMENSION(ncols,nrows),INTENT(OUT) :: gam_l

      gam_l = (0.49*lai) / ( (1+0.2*(lai**2))**0.5 )

      RETURN
      END SUBROUTINE GAMMA_LAI
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!.....2) Calculate GAM_P (GAMMA_P)
!-----------------------------------------------------------------------
!             GAMMA_P = 0.0         a<=0, a>=180, sin(a) <= 0.0
!
!             GAMMA_P = sin(a)[ 2.46*(1+0.0005(Pdaily-400))*PHI - 0.9*PHI^2 ]
!                                   0<a<180, sin(a) > 0.0
!           where PHI    = above canopy PPFD transmission (non-dimension)
!                 Pdaily = daily average above canopy PPFD (umol/m2s)
!                 a      = solar angle (degree)
!
!                 Note: AAA = 2.46*BBB*PHI - 0.9*PHI^2
!                       BBB = (1+0.0005(Pdaily-400))
!                       GAMMA_P = sin(a)*AAA
!
!                       Pac
!             PHI = -----------
!                   sin(a)*Ptoa
!           where Pac  = above canopy PPFD (umol/m2s)
!                 Ptoa = PPFD at the top of atmosphere (umol/m2s)
!
!             Pac =  SRAD * 4.766 mmmol/m2-s * 0.5
!
!             Ptoa = 3000 + 99*cos[2*3.14-( DOY-10)/365 )]
!           where DOY = day of year
!
!     SUBROUTINE GAMMA_P returns the GAMMA_P values
!-----------------------------------------------------------------------
      SUBROUTINE GAMMA_P( JDATE, JTIME, LAT, LONG,          
     &                    NCOLS, NROWS, PPFD, D_PPFD, GAM_P )

      IMPLICIT NONE

! input
      INTEGER,INTENT(IN) :: NCOLS, NROWS, JDATE, JTIME

      REAL,DIMENSION(NCOLS,NROWS),INTENT(IN) :: LAT, LONG   
      ! Photosynthetic Photon Flux Density: instantaneous, daily
      REAL,DIMENSION(NCOLS,NROWS),INTENT(IN) :: PPFD, D_PPFD   
! output
      REAL,DIMENSION(NCOLS,NROWS),INTENT(OUT) :: GAM_P          ! GAMMA_P

! Local parameters
      CHARACTER*256 MESG          ! message buffer
      CHARACTER*16  :: FUNCNAME = 'GAMMA_P'

      INTEGER,DIMENSION(NCOLS,NROWS) :: DAY    ! DAY is DOY (JDATE)
      REAL ::   AAA, BBB
      REAL ::   BETA                           ! Solar zenith angle
      REAL ::   SINbeta                        ! sin(beta)
      REAL,DIMENSION(NCOLS,NROWS) :: HOUR      ! HOUR is solar hour
      REAL,DIMENSION(NCOLS,NROWS) :: SZANGLE   ! Solar zenith angle array
      REAL,DIMENSION(NCOLS,NROWS) :: Ptoa, Pac, PHI

      INTEGER :: I, J

!...  Constants
      REAL,PARAMETER :: PI = 3.14159, D2RAD = PI/180.0, RAD2D = 180.0/PI

!...  Convert date and time format to local time
      ! DAY is Julian day
      DAY(:,:)  =  MOD(JDATE,1000)

      ! Convert from XXXXXX format to XX.XX (solar hour)
      ! HOUR = 0 -> 23.xx
      ! Solar hour
      HOUR = JTIME/10000 + LONG/15
            
      WHERE ( HOUR .LT. 0.0 )
        HOUR = HOUR + 24.0
        DAY = DAY - 1
      ENDWHERE

!...  Begin estimating gamma_p
      ! getting solar radiation
      Pac = PPFD

!...  Initialize parameters
      SZANGLE = 0.
      Ptoa = 0.
      PHI = 0.

      DO I = 1, NCOLS
        DO J = 1, NROWS
          ! Get solar elevation angle
          CALL SOLARANGLE(DAY(I,J),HOUR(I,J),LAT(I,J),SINbeta)
          BETA = ASIN(SINbeta)*RAD2D            ! Degree
          SZANGLE(I,J) =  BETA                  ! Degree
          IF (SINbeta .LE. 0.0) THEN
            GAM_P(I,J) = 0.0
          ELSEIF (SINbeta .GT. 0.0) THEN
            Ptoa(I,J) = 3000.0 + 99.0 *COS(2*3.14*(DAY(I,J)-10)/365)

            PHI(I,J) = Pac(I,J)/(SINbeta * Ptoa(I,J))

            BBB = 1 + 0.0005*( D_PPFD(I,J)-400 )
            AAA = ( 2.46 * BBB * PHI(I,J) ) - ( 0.9 * PHI(I,J)**2 )

            GAM_P(I,J) = SINbeta * AAA
          ELSE
            MESG = 'Error: Solar angle is invalid'
            CALL M3EXIT(FUNCNAME,JDATE,JTIME,MESG,2)

          ENDIF

          ! Screening the unforced errors
          ! IF solar elevation angle is less than 1 THEN
          ! gamma_p can not be greater than 0.1.
          IF (BETA .LT. 1.0 .AND. GAM_P(I,J) .GT. 0.1) THEN
            GAM_P(I,J) = 0.0
          ENDIF

        ENDDO    ! End loop for NROWS
      ENDDO      ! End loop for NCOLS

      RETURN
      END SUBROUTINE GAMMA_P
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!.....3) Calculate GAM_T (GAMMA_T) for isoprene
!-----------------------------------------------------------------------
!                          Eopt*CT2*exp(CT1*x)
!             GAMMA_T =  ------------------------
!                        [CT2-CT1*(1-exp(CT2*x))]
!           where x      = [ (1/Topt)-(1/Thr) ] / 0.00831
!                 Eopt   = 1.75*exp(0.08(Tdaily-297)
!                 CT1    = 80
!                 CT2    = 200
!                 Thr    = hourly average air temperature (K)
!                 Tdaily = daily average air temperature (K)
!                 Topt   = 313 + 0.6(Tdaily-297)
!
!                 Note: AAA = Eopt*CT2*exp(CT1*x)
!                       BBB = [CT2-CT1*(1-exp(CT2*x))]
!                       GAMMA_T = AAA/BBB
!
!     SUBROUTINE GAMMA_TLD returns the GAMMA_T value for isoprene
!-----------------------------------------------------------------------
      SUBROUTINE GAMMA_TLD( NCOLS, NROWS, TEMP, 
     &             D_TEMP, GAM_T,SPC_NAME )

      IMPLICIT NONE
      INCLUDE 'EACO.EXT'
      INTEGER,EXTERNAL :: INDEX1
! input
      INTEGER,INTENT(IN) :: ncols, nrows
      REAL,DIMENSION(ncols,nrows),INTENT(IN) :: TEMP,D_TEMP   ! daily, hourly surface temperature 
! output
      REAL,DIMENSION(ncols,nrows),INTENT(OUT) :: GAM_T       ! GAMMA_T
       CHARACTER*16,INTENT(IN) :: SPC_NAME
       REAL,PARAMETER :: CT2 =200.0  
! Local parameters
        INTEGER :: SPCNUM
      REAL,DIMENSION(ncols,nrows) :: Eopt, Topt, X, AAA, BBB
   
       SPCNUM = INDEX1(SPC_NAME,N_MGN_SPC,MGN_SPC)

      Eopt = Cceo(SPCNUM) * exp(0.08*(D_TEMP-297.0))
      Topt = 313.0 + ( 0.6*(D_TEMP-297.0) )
      X = ( (1/Topt)-(1/TEMP) ) / 0.00831

      AAA = Eopt*CT2*exp(CT1(SPCNUM)*X)
      BBB = ( CT2-CT1(SPCNUM)*( 1-exp(CT2*X) ) )
      GAM_T = AAA/BBB

      RETURN
      END SUBROUTINE GAMMA_TLD
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!.....4) Calculate GAM_T (GAMMA_T) for non-isoprene
!-----------------------------------------------------------------------
!
!             GAMMA_T =  exp[TDP_FCT*(T-Ts)]
!           where TDP_FCT = temperature dependent parameter ('beta')
!                 Ts     = standard temperature (normally 303K, 30C)
!
!     SUBROUTINE GAMMA_TLI returns the GAMMA_T value for non-isoprene
!-----------------------------------------------------------------------
      SUBROUTINE GAMMA_TLI( NCOLS, NROWS, SPCNAM, TEMP, GAM_T)

      IMPLICIT NONE
     
      INCLUDE 'EACO.EXT'

      INTEGER,EXTERNAL :: INDEX1
      CHARACTER*16  SPCNAM
      INTEGER :: NCOLS, NROWS
      INTEGER :: SPCNUM                             ! Species number
      REAL,DIMENSION(NCOLS,NROWS) :: TEMP, GAM_T
      REAL, PARAMETER :: Ts = 303.0

!--end of declarations--

      SPCNUM = INDEX1(SPCNAM,N_TDF_SPC,TDF_SPC)
      GAM_T = exp( TDF_PRM(SPCNUM)*(TEMP-Ts) )

      RETURN
      END SUBROUTINE GAMMA_TLI
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!.....5) Calculate GAM_A (GAMMA_age)
!-----------------------------------------------------------------------
!
!             GAMMA_age = Fnew*Anew + Fgro*Agro + Fmat*Amat + Fold*Aold
!           where Fnew = new foliage fraction
!                 Fgro = growing foliage fraction
!                 Fmat = mature foliage fraction
!                 Fold = old foliage fraction
!                 Anew = relative emission activity for new foliage
!                 Agro = relative emission activity for growing foliage
!                 Amat = relative emission activity for mature foliage
!                 Aold = relative emission activity for old foliage
!
!
!             For foliage fraction
!             Case 1) LAIc = LAIp
!             Fnew = 0.0  , Fgro = 0.1  , Fmat = 0.8  , Fold = 0.1
!
!             Case 2) LAIp > LAIc
!             Fnew = 0.0  , Fgro = 0.0
!             Fmat = 1-Fold
!             Fold = (LAIp-LAIc)/LAIp
!
!             Case 3) LAIp < LAIc
!             Fnew = 1-(LAIp/LAIc)                       t <= ti
!                  = (ti/t) * ( 1-(LAIp/LAIc) )          t >  ti
!
!             Fmat = LAIp/LAIc                           t <= tm
!                  = (LAIp/LAIc) +
!                      ( (t-tm)/t ) * ( 1-(LAIp/LAIc) )  t >  tm
!
!             Fgro = 1 - Fnew - Fmat
!             Fold = 0.0
!
!           where
!             ti = 5 + (0.7*(300-Tt))                   Tt <= 303
!                = 2.9                                  Tt >  303
!             tm = 2.3*ti
!
!             t  = length of the time step (days)
!             ti = number of days between budbreak and the induction of
!                  emission
!             tm = number of days between budbreak and the initiation of
!                  peak emissions rates
!             Tt = average temperature (K) near top of the canopy during
!                  current time period (daily ave temp for this case)
!
!
!             For relative emission activity
!             Case 1) Constant
!             Anew = 1.0  , Agro = 1.0  , Amat = 1.0  , Aold = 1.0
!
!             Case 2) Monoterpenes
!             Anew = 2.0  , Agro = 1.8  , Amat = 0.95 , Aold = 1.0
!
!             Case 3) Sesquiterpenes
!             Anew = 0.4  , Agro = 0.6  , Amat = 1.075, Aold = 1.0
!
!             Case 4) Methanol
!             Anew = 3.0  , Agro = 2.6  , Amat = 0.85 , Aold = 1.0
!
!             Case 5) Isoprene
!             Anew = 0.05 , Agro = 0.6  , Amat = 1.125, Aold = 1.0
!
!     SUBROUTINE GAMMA_A returns GAMMA_A
!-----------------------------------------------------------------------
      SUBROUTINE GAMMA_A( JDATE, JTIME, NCOLS, NROWS, SPC_NAME,  
     &                    LAIARp, LAIARc, TSTLEN, D_TEMP, GAM_A )

      IMPLICIT NONE

      INCLUDE 'EACO.EXT'

      INTEGER,EXTERNAL :: INDEX1
! input
      INTEGER,INTENT(IN) :: NCOLS, NROWS, JDATE, JTIME,TSTLEN
      CHARACTER*16,INTENT(IN) :: SPC_NAME
      REAL,DIMENSION(ncols,nrows),INTENT(IN) :: D_TEMP, LAIARp, LAIARc
! output
      REAL,DIMENSION(ncols,nrows),INTENT(OUT) :: GAM_A

! Local parameters
      CHARACTER*16  :: FUNCNAME = 'GAMMA_A'

      INTEGER :: AINDX                            ! relative emission acitivity index
      INTEGER :: SPCNUM
      INTEGER :: t                                ! time step
      CHARACTER*256 :: MESG                       ! message buffer

      REAL,DIMENSION(ncols,nrows) :: Fnew, Fgro, Fmat, Fold
      REAL,DIMENSION(ncols,nrows) :: LAIp,LAIc    ! LAI at previous,current time steps
      REAL,DIMENSION(ncols,nrows) :: ti,tm        ! number of days between budbreak
                                                  ! and induction of emission, 
                                                  ! initiation of peak emissions rates
      REAL,DIMENSION(ncols,nrows) :: Tt           ! daily average temperature (K)

!...  Choose relative emission activity
!--------code by Xuemei Wang 11/04/2007----------------       
      SPCNUM = INDEX1(SPC_NAME,N_MGN_SPC,MGN_SPC)
      AINDX =  REA_INDEX(SPCNUM)
!---------------------------------------------------        
! local parameter arrays
      t = TSTLEN
      LAIc = LAIARc
      LAIp = LAIARp
      Tt   = D_TEMP

!... Calculate foliage fraction

      WHERE (LAIp .LT. LAIc) 

!        Calculate ti and tm
        WHERE (Tt .LE. 303.0) 
          ti = 5.0 + 0.7*(300-Tt)
        ELSEWHERE (Tt .GT. 303.0) 
          ti = 2.9
        ENDWHERE
        tm = 2.3*ti
 
!       Calculate Fnew and Fmat, then Fgro and Fold
!       Fnew
        WHERE (ti .GE. t) 
          Fnew = 1.0 - (LAIp/LAIc)
        ELSEWHERE (ti .LT. t) 
          Fnew = (ti/t) * ( 1-(LAIp/LAIc) )
        ENDWHERE
 
!       Fmat
        WHERE (tm .GE. t) 
          Fmat = LAIp/LAIc
        ELSEWHERE (tm .LT. t)
          Fmat = (LAIp/LAIc) + ( (t-tm)/t ) * ( 1-(LAIp/LAIc) )
        ENDWHERE

        Fgro = 1.0 - Fnew - Fmat
        Fold = 0.0
         
      ELSEWHERE (LAIp .EQ. LAIc) 

        Fnew = 0.0
        Fgro = 0.1
        Fmat = 0.8
        Fold = 0.1

      ELSEWHERE (LAIp .GT. LAIc) 
 
        Fnew = 0.0
        Fgro = 0.0
        Fold = ( LAIp-LAIc ) / LAIp
        Fmat = 1-Fold

      ENDWHERE

!...  Calculate GAMMA_A
      GAM_A = Fnew*Anew(AINDX) + Fgro*Agro(AINDX) +  
     &         Fmat*Amat(AINDX) + Fold*Aold(AINDX)

      RETURN
      END SUBROUTINE GAMMA_A

!-----------------------------------------------------------------------
!.....6) Calculate GAM_SMT (GAMMA_SM)
!-----------------------------------------------------------------------
!
!             GAMMA_SM =     1.0   (non-dimension)
!
!
!     SUBROUTINE GAMMA_S returns the GAMMA_SM values
!-----------------------------------------------------------------------
      SUBROUTINE GAMMA_S( NCOLS, NROWS, GAM_S )

      IMPLICIT NONE

      INTEGER :: NCOLS, NROWS
      REAL,DIMENSION(NCOLS,NROWS) :: GAM_S

      GAM_S = 1.0

      RETURN
      END SUBROUTINE GAMMA_S

!=======================================================================
!-----------------------------------------------------------------------
!.....7) Calculate GAM_CO2(GAMMA_CO2)
!-----------------------------------------------------------------------
!
!             GAMMA_CO2 =     1.0   (non-dimension)
!             When CO2 =400ppm
!
!     SUBROUTINE GAM_CO2 returns the GAMMA_CO2 values
!    Xuemei Wang-2009-06-22 
!-----------------------------------------------------------------------
      SUBROUTINE GAMMA_CO2( NCOLS, NROWS,CO2, GAM_CO2 )

      IMPLICIT NONE

      INTEGER :: NCOLS, NROWS
      REAL,DIMENSION(NCOLS,NROWS) :: GAM_CO2,CO2,Ci
      REAL, PARAMETER :: ISmax = 1.344, h=1.4614
      REAL, PARAMETER :: Cstar =585
    
   
         Ci = 0.7* CO2
         WHERE (CO2.eq.400.0)
       GAM_CO2 = 1.0
       ELSEWHERE
       GAM_CO2 = ISmax- ((ISmax*Ci**h) /(Cstar**h+Ci**h))
         ENDWHERE

      RETURN
      END SUBROUTINE GAMMA_CO2
     
!=======================================================================
!=======================================================================
!-----------------------------------------------------------------------
!.....8) Calculate GAMMA_LAIbidir(gam_LAIbidir,LAI)
!-----------------------------------------------------------------------
!From Alex Guenther 2010-01-26
!If lai < 2 Then
!gammaLAIbidir= 0.5 * lai
!ElseIf lai <= 6 Then
!gammaLAIbidir= 1 - 0.0625 * (lai - 2)
!Else
!gammaLAIbidir= 0.75
!End If
!
!     SUBROUTINE GAMMA_LAIbidir returns the gam_LAIbidir values
!    Xuemei Wang-2010-01-28
!
!-----------------------------------------------------------------------
       SUBROUTINE GAMMA_LAIbidir(NCOLS, NROWS,LAI,GAM_LAIbidir)

          IMPLICIT NONE

       INTEGER,INTENT(IN) :: NCOLS, NROWS
       INTEGER :: I,J
       REAL,DIMENSION(NCOLS, NROWS),INTENT(IN) ::  LAI
       REAL,DIMENSION(NCOLS, NROWS),INTENT(OUT) :: GAM_LAIbidir
        DO I = 1,NCOLS
        DO J = 1, NROWS
         IF(LAI(I,J) < 2) THEN
        GAM_LAIbidir =  0.5 * LAI
        ELSEIF (LAI(I,J) .LE. 6 .AND. LAI(I,J) .GE. 2) THEN 
        GAM_LAIbidir = 1 - 0.0625 * (LAI(I,J) - 2)
        ELSE
        GAM_LAIbidir = 0.75
        ENDIF
     
        ENDDO
        ENDDO

        RETURN
      END  SUBROUTINE GAMMA_LAIbidir
!=======================================================================

      END MODULE GAMMA_etc
