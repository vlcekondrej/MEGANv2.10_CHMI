!   Modified by Xuemei Wang--Nov. 2007
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!   Input and output files must be selected before starting the program
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!!
!   Input varibles
!
!   Day                  Julian day
!   Lat                  Latitude
!   Hour                 Hour of the day
!   Tc                   Temperature [C]
!   PPFD           Incoming photosynthetic active radiation [umol/m2/s1]
!   Wind                 Wind speed [m s-1]
!   Humidity             Relative humidity [%]
!   Cantypye             Defines set of canopy characteristics
!   LAI                  Leaf area index [m2 per m2 ground area]
!   DI                   ???
!   Pres                 Pressure [Pa]
!
!   Used variables:
!
!   PPFDfrac             Fraction of total solar radiation that is PPFD
!   Solar                Solar radiation [W/m2]
!   Maxsolar             Maximum of solar radiation
!   Beta                 Sin of solar angle above horizon
!   Sinbeta              Solar angle above horizon
!   TairK0               Above canopy air temperature [K]
!       TairK            Array of canopy air temperature [K]
!   Ws0                  Above canopy wind speed [m/s]
!   Ws                   Array of canopy wind speed [m/s]
!   HumidairPA0          Above canopy ambient humidity [Pa]
!       HumidairPa       Array of canopy ambient humidity in [Pa]
!   StomataDI            Index for water status of leaves. used to modify stomatal conductance
!   Transmis             Transmission of PPFD that is diffuse
!   Difffrac             Fraction of PPFD that is diffuse
!   PPFDfrac             Fraction of solar rad that is PPFD
!   Trate                Stability of boundary ???
!   SH                   Sensible heat flux ???
!       VPgausWt         Array of gaussian weighting factors
!   VPgausDis            Array of gaussian weighting factors
!   VPslwWT              Array of gaussian weighting factors
!   SunFrac              Array of the fraction of sun leaves. i = 1 is the top canopy layer, 2 is the next layer, etc.
!   SunPPFD              Array of incoming (NOT absorbed) PPFD on a sun leaf [umol/m2/s]
!   ShadePPFD            Array of incoming (NOT absorbed) PPFD on a shade leaf [umol/m2/s]
!   SunQv                Array of visible radiation (in and out) fluxes on sun leaves
!   ShadeQv              Array of absorbed visible radiation (in and out) fluxes on shade leaves
!   SunQn                Array of absorbed near IR radiation (in and out) fluxes on sun leaves
!   ShadeQn              Array of absorbed near IR radiation (in and out) fluxes on shade leaves
!   SunleafTK            Array of leaf temperature for sun leaves [K]
!   SunleafSH            Array of sensible heat flux for sun leaves [W/m2]
!   SunleafLH            Array of latent heat flux for sun leaves [W/m2]
!   SunleafIR            Array of infrared flux for sun leaves [W/m2]
!   ShadeleafTK          Array of leaf temperature for shade leaves [K]
!   ShadeleafSH          Array of sensible heat flux for shade leaves [W/m2]
!   ShadeleafLH          Array of latent heat flux for shade leaves [W/m2]
!   ShadeleafIR          Array of infrared flux for shade leaves [W/m2]
!   QbAbsV, QbAbsN       Absorbed direct beam light for visible and near infra red
!   QdAbsV, QdAbsN       Array of absorbed diffuse light for visible and near infra red
!   QsAbsV, QsAbsN       Array of absorbed scattered light for visible and near infra red
!   QBeamV, QBeamN       Above canopy beam (direct) light for visible and near infra red
!   QDiffV, QDiffN       Above canopy diffuse light for visible and near infra red
!       Ea1pLayer        Array of emission activity of light per layer
!       Ea1tLayer        Array of emission activity of temperature per layer
!       Ea1Layer         Array of companied emission activity
!       Ea1pCanopy       Total emission activity of light
!       EatiLayer        Array of emission activity of temperature indendent per layer
!   Ea1tCanopy           Total emission activity of temperature depedent factor
!   Ea1Canopy            Total companied emission activity
!   EatiCanopy           Total emission activity of temperature indepedent factor
!   Calcbeta             Function: Calculation of solar zenith angle
!   WaterVapPres         Function: Convert water mixing ratio (kg/kg) to water vapor pressure
!   Stability            Function: Temperature lapse rate
!   Ea1t99               Function: Temperature dependence activity factor for emission type 1
!   Ea1p99               Function: Light dependence activity factor for emission
!   Ealti                Function: Temperature independence activity factor for emission
!   DIstomata            Function:
!   CalcEccentricity     Function:
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      SUBROUTINE GAMME_CE(JDATE,JTIME,LAT, LONG,Tc,T24,T240,          
     &                     PPFD0, PPFD24, PPFD240, Wind,              
     &             Humidity,Cantype,LAI,Pres,DI,NrCha,NrTyp,         
     &             Canopychar,SPCNAME,Ea1Canopy,EatiCanopy)
     
 
      IMPLICIT NONE
      INCLUDE 'CONST_CANOPY.EXT'

! input
      INTEGER,INTENT(IN) :: JDATE,JTIME,NrCha,NrTyp
      REAL,INTENT(IN) :: LONG,LAT,Cantype
      REAL,INTENT(IN) :: Tc, Pres, Wind,Humidity, LAI, DI,            
     &                    PPFD0,PPFD24,PPFD240,T24,T240
      ! Array of canopy characteristics for NrTyp of canopy type
      REAL,DIMENSION(NrCha,NrTyp),INTENT(IN) :: Canopychar             
      CHARACTER*16,INTENT(IN) :: SPCNAME

! output
      REAL,INTENT(OUT) ::  Ea1Canopy,EatiCanopy 

! local variables
      INTEGER :: i, JJ, KK, K, Day

      REAL :: Sinbeta,Beta
      REAL,DIMENSION(Layers) ::  VPgausWt, VPgausDis2, 
     &  VPgausDis, VPslwWT, Sunfrac, QdAbsV, QsAbsV, QdAbsn,     
     &  QsAbsn, SunQv, ShadeQv, SunQn, ShadeQn, SunPPFD,                
     &  ShadePPFD, TairK, HumidairPa, Ws, Sunleaftk, SunleafSH,            
     &  SunleafLH,SunleafIR, Shadeleaftk, ShadeleafSH,            
     &  ShadeleafLH,ShadeleafIR,Ea1pLayer, Ea1tLayer, Ea1Layer,  
     &  EatiLayer

      REAL :: Hour, Solar, Maxsolar,  
     &         Difffrac,PPFD, PPFDfrac, QbAbsn,                  
     &         Trate, StomataDI, Qbeamv,Qdiffv, Qbeamn, Qdiffn,  
     &         QbAbsV,Ea1tCanopy, Ea1pCanopy,                    
     &         TairK0, HumidairPa0,Ws0, SH
                       
      REAL ::  DIstomata, CalcEccentricity,WaterVapPres,        
     &          Stability, Ea1t99, Ea1p99,Calcbeta,Ealti99

!---------------------------Header over--------------------------------

      Day  =  MOD(JDATE,1000)
! Convert from XXXXXX format to XX.XX (solar hour)
      ! HOUR = 0 -> 23.xx
! Solar hour
      Hour  = JTIME/10000 + LONG /15
!        Hour  = JTIME + LONG /15
      IF ( Hour  .LT. 0.0 ) THEN
        Hour  = Hour  + 24.0
        Day  = Day  - 1
      ELSEIF ( Hour  .GT. 24.0 ) THEN
        print*,'Invalid hour: HOUR  is ', Hour 
      ENDIF

      Beta   = Calcbeta(Day , Lat , Hour )
      Sinbeta    = SIN(Beta  / 57.29578)
      TairK0     = Tc   !See the input unit(K)
      Ws0        = Wind 
      PPFD  = PPFD0 
      StomataDI  = DIstomata(DI )
!      Solar      = PPFD /ConvertWm2toUmolm2s*2
!      Solar      = PPFD/2.1
       Solar      = PPFD/2.25
      Maxsolar   = Sinbeta  * SolarConstant * CalcEccentricity(Day )

!	print*,Solar,Hour,Sinbeta 
      Call GaussianIntegration(VPgausWt, VPgausDis,VPgausDis2, Layers)
   
      Call SolarFractions(Solar, Maxsolar, Qdiffv,Qbeamv,Qdiffn,Qbeamn)

!      Qdiffv = PPFDfrac * Solar * Difffrac
!      Qbeamv = PPFDfrac * Solar * (1 - Difffrac)
!      Qdiffn = (1 - PPFDfrac ) * Solar * Difffrac 
!See the input unit(K) 
!      Qbeamn = (1 - PPFDfrac ) * Solar * (1 - Difffrac )
!	 print*,Qdiffv,Qbeamv
      Call WeightSLW(VPgausDis, VPgausWt, LAI, Layers, VPslwWT)
         
      Call CanopyRad(VPgausDis, Layers, LAI, Sinbeta, Qbeamv,    
     &     Qdiffv, Qbeamn, Qdiffn,Cantype ,Canopychar, Sunfrac,  
     &     QbAbsV, QdAbsV, QsAbsV, QbAbsn, QdAbsn, QsAbsn, SunQv, 
     &     ShadeQv, SunQn, ShadeQn, SunPPFD, ShadePPFD, NrCha,    
     &              NrTyp)


      HumidairPa0  =  WaterVapPres(Humidity , Pres , WaterAirRatio)
      Trate    =  Stability(Canopychar, Cantype, Solar , NrCha, NrTyp)

      Call CanopyEB(Trate, Layers, VPgausDis, Canopychar, Cantype,      
     &                StomataDI, TairK, HumidairPa, Ws, SunPPFD,         
     &                ShadePPFD, SunQv, ShadeQv, SunQn, ShadeQn,         
     &                Sunleaftk, SunleafSH, SunleafLH, SunleafIR,        
     &                Shadeleaftk,ShadeleafSH,ShadeleafLH,ShadeleafIR,   
     &                NrCha, NrTyp, Ws0, TairK0, HumidairPa0)
!        print*, Trate, Shadeleaftk(1), Sunleaftk(1)

      Ea1tCanopy = 0.
      Ea1pCanopy = 0.
      Ea1Canopy  = 0.
      EatiCanopy = 0.
      SH         = 0.
      
      DO i=1,layers    

        Ea1tLayer(i) = Ea1t99(Sunleaftk(i), T24, T240,SPCNAME) *
     &                  Sunfrac(i) + 
     &                 Ea1t99(Shadeleaftk(i),T24,T240,SPCNAME) *
     &                   (1-Sunfrac(i))

! pstd = 200 for sun leaves
! pstd = 50 for shade leaves
        Ea1pLayer(i) = Ea1p99(SunPPFD(i),PPFD24*0.5,           
     &                        PPFD240*0.5,pstd_Sun)* Sunfrac(i) +       
     &                 Ea1p99(ShadePPFD(i),PPFD24*0.16,        
     &                     PPFD240*0.16,pstd_Shade)*(1-Sunfrac(i))  

        Ea1Layer(i)  = Ea1t99(Sunleaftk(i), T24, T240,SPCNAME) *       
     &                 Ea1p99(SunPPFD(i),PPFD24*0.5,           
     &                        PPFD240*0.5,pstd_Sun) * Sunfrac(i) +   
     &                 Ea1t99(Shadeleaftk(i), T24, T240,SPCNAME) *     
     &                 Ea1p99(ShadePPFD(i),PPFD24*0.16,        
     &                     PPFD240*0.16,pstd_Shade)*(1-Sunfrac(i))
          
        EatiLayer(i) = Ealti99(SPCNAME,Sunleaftk(i))* Sunfrac(i)+   
     &              Ealti99(SPCNAME,Shadeleaftk(i))*(1-Sunfrac(i))
 
    
      ENDDO

      Ea1pCanopy = SUM( Ea1pLayer(:) * VPslwWT(:) * VPgausWt(:) )
      Ea1tCanopy = SUM( Ea1tLayer(:) * VPslwWT(:) * VPgausWt(:) )
      Ea1Canopy  = SUM( Ea1Layer(:)  * VPslwWT(:) * VPgausWt(:) )
      EatiCanopy = SUM( EatiLayer(:) * VPslwWT(:) * VPgausWt(:) )
      Ea1Canopy = Ea1Canopy *Cce* LAI

! this quantity is apparently not passed out of the subroutine
!      SH = SUM( (SunleafSH(:) * Sunfrac(:) +                 
!     &           ShadeleafSH(:) * (1 - Sunfrac(:))) *       
!     &           LAI  * VPgausWt(:)                   )

       RETURN
       END SUBROUTINE  GAMME_CE

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   SUBROUTINE GaussianIntegration
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      SUBROUTINE GaussianIntegration
     &           (Weightgauss, Distgauss,Distgauss2, Layers)
 
      IMPLICIT NONE
      INTEGER ,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)

      INTEGER,INTENT(IN) ::  Layers

      REAL,DIMENSION(Layers),INTENT(OUT) :: 
     &                         Weightgauss, Distgauss,Distgauss2

! local variables
      INTEGER ::  i
!--------------------------------------------------------------------

      IF (Layers .EQ. 1) THEN
        Weightgauss(1) = 1
        Distgauss(1)   = 0.5
        Distgauss2(1)   = 1 
      ELSEIF (Layers .EQ. 3) THEN
        Weightgauss(1) = 0.277778
        Weightgauss(2) = 0.444444
        Weightgauss(3) = 0.277778
        Distgauss(1)   = 0.112702
        Distgauss(2)   = 0.5
        Distgauss(3)   = 0.887298
        Distgauss2(1)   = 0.277778
        Distgauss2(2)   = 0.722222
        Distgauss2(3)   = 1
      ELSEIF (Layers .EQ. 5) THEN
        Weightgauss(1) = 0.1184635
        Weightgauss(2) = 0.2393144
        Weightgauss(3) = 0.284444444
        Weightgauss(4) = 0.2393144
        Weightgauss(5) = 0.1184635
        Distgauss(1)   = 0.0469101
        Distgauss(2)   = 0.2307534
        Distgauss(3)   = 0.5
        Distgauss(4)   = 0.7692465
        Distgauss(5)   = 0.9530899
        Distgauss2(1)   = 0.1184635
        Distgauss2(2)   = 0.3577778
        Distgauss2(3)   = 0.6422222
        Distgauss2(4)   = 0.881536
        Distgauss2(5)   = 1.0

      ELSE
        DO i = 1, Layers
          Weightgauss(i) = 1. / Layers
          Distgauss(i)   = (i - 0.5) / Layers
         Distgauss2(i) = i/Layers
        ENDDO
      ENDIF

      RETURN
      END SUBROUTINE GaussianIntegration

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   SUBROUTINE WeightSLW
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      SUBROUTINE WeightSLW(Distgauss, Weightgauss, LAI, Layers, SLW)

      IMPLICIT NONE
      INTEGER,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
  
      INTEGER,INTENT(IN) ::  Layers
      REAL,INTENT(IN) ::  LAI
      REAL,DIMENSION(Layers),INTENT(IN) :: Distgauss, Weightgauss

      REAL,DIMENSION(Layers),INTENT(OUT) :: SLW

! local variables
      INTEGER :: i
!--------------------------------------------------           

      DO i = 1, Layers
       SLW(i) = 0.63 + 0.37 * EXP(-((LAI  *
     &            Distgauss(i)) - 1))

      ENDDO


      RETURN
      END SUBROUTINE WeightSLW

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   SUBROUTINE SolarFractions
!   Transmission, fraction of PPFD that is diffuse, 
!        fraction of solar rad that is PPFD
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      SUBROUTINE SolarFractions( Solar, Maxsolar,  
     &                          Qdiffv,Qbeamv, Qdiffn, Qbeamn)

      IMPLICIT NONE
      INTEGER,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
 
!      INTEGER,INTENT(IN) :: Timeperiod
      REAL,INTENT(IN) :: Solar, Maxsolar

      REAL,INTENT(OUT) ::  Qdiffv,Qbeamv, Qdiffn, Qbeamn
      REAL :: FracDiff, PPFDfrac,PPFDdifFrac,Qv, Qn

! internal variables
      INTEGER :: I,J
      REAL ::  Transmis 
!-----------------------------------------------------

!      IF (Timeperiod .EQ. 1) THEN     ! Daily transmission
!        TransMin   = 0.26
!        TransSlope= 1.655
!      ELSE                       ! Hourly transmission
!        TransMin   = 0.26
!        TransSlope = 1.655
!      ENDIF
         
      IF (Maxsolar  <= 0) THEN
        Transmis  = 0.5
       ELSEIF (Maxsolar < Solar) THEN
      Transmis  = 1.0
        ELSE
       Transmis  = Solar  / Maxsolar 
      ENDIF

! Estimate diffuse fraction based on daily transmission (Roderick 1999, Goudriann and Van Laar 1994- P.33)

!      IF (Transmis  > 0.81) THEN
!        FracDiff  = 0.05
!      ELSEIF (Transmis  > TransMin) THEN
!        FracDiff  = 0.96-TransSlope * (Transmis  - TransMin)
!      ELSE
!        FracDiff  = 0.96
!      ENDIF

! The fraction of total solar radiation that is PPFD (43% to 55%) 
! G. and L. 84
!      PPFDfrac  = 0.43 + FracDiff  * 0.12

!FracDiff is based on Lizaso 2005
!modified by Xuemei 2010-01-26 according to Alex's document
        FracDiff = 0.156 + 0.86/(1 + EXP(11.1*(Transmis -0.53)))

!PPFDfrac is based on G.L. 84
!modified by Xuemei 2010-01-26 according to Alex's document
        PPFDfrac  = 0.55 -Transmis*0.12

!PPFDdifFrac is based on data in Jacovides 2007
!modified by Xuemei 2010-01-26 according to Alex's document
        PPFDdifFrac = FracDiff *(1.06 + Transmis*0.4)  

! Calculte  Qdiffv,Qbeamv, Qdiffn, Qbeamn in the subroutine
! modified by Xuemei 2010-01-26 according to Alex's document
        IF (PPFDdifFrac > 1.0) THEN
        PPFDdifFrac  = 1.0
        ENDIF
 
        Qv  = PPFDfrac * Solar
        Qdiffv = Qv * PPFDdifFrac
        Qbeamv = Qv - Qdiffv
        Qn = Solar - Qv
        Qdiffn =  Qn * FracDiff
        Qbeamn =  Qn - Qdiffn
            

      RETURN
      END SUBROUTINE SolarFractions

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   Subroutine CanopyRad
!
!   Canopy light environment model
!   Code developed by Alex Guenther, based on Spitters et al. (1986), 
!   Goudrian and Laar (1994), Leuning (1997)
!   Initial code 8-99, modified 7-2000 and 12-2001
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      SUBROUTINE CanopyRad(Distgauss, Layers, LAI, Sinbeta,   
     &            Qbeamv, Qdiffv, Qbeamn, Qdiffn, Cantype,     
     &            Canopychar, Sunfrac, QbAbsV, QdAbsV, QsAbsV, 
     &            QbAbsn, QdAbsn, QsAbsn, SunQv,               
     &            ShadeQv, SunQn, ShadeQn, SunPPFD, ShadePPFD, 
     &            NrCha, NrTyp)

      IMPLICIT NONE

! input
      INTEGER,INTENT(IN) :: Layers, NrCha, NrTyp, Cantype
      REAL,INTENT(IN) :: Qbeamv,Qdiffv,Sinbeta,LAI,Qbeamn,Qdiffn
      REAL,DIMENSION(Layers),INTENT(IN) :: Distgauss

! output
      REAL,INTENT(OUT) :: QbAbsV, QbAbsn

      REAL,DIMENSION(Layers),INTENT(OUT) :: ShadePPFD, SunPPFD,  
     &                QdAbsv, QsAbsv, QsAbsn, ShadeQv,  SunQn,    
     &                QdAbsn, SunQv, ShadeQn, Sunfrac

      REAL,DIMENSION(NrCha,NrTyp),INTENT(OUT) :: Canopychar

! internal variables
      INTEGER :: i
      REAL :: ScatV, ScatN, RefldV, RefldN, ReflbV, ReflbN,      
     &         Kb, Kd, KbpV, KbpN, KdpV, KdpN, LAIdepth, Cluster, 
     &         QdAbsVL, QsAbsVL, QdAbsNL, QsAbsNL

! Stefan-boltzman constant  W m-2 K-4
      REAL,PARAMETER :: Sb = 0.0000000567    
!     REAL,PARAMETER :: ConvertPPFD = 4.766    
      REAL,PARAMETER :: ConvertShadePPFD = 4.6
      REAL,PARAMETER :: ConvertSunPPFD = 4.0
!---------------------------------------------------------------------
        
      IF (((Qbeamv  + Qdiffv ) > 0.001) .AND.    
     &     (Sinbeta  > 0.00002) .AND.             
     &     (LAI  > 0.001)) THEN       ! Daytime

! Scattering coefficients (scatV,scatN), diffuse and beam reflection 
! coefficients (ref..) for visible or near IR
        ScatV   = Canopychar(5,Cantype)
        ScatN   = Canopychar(6,Cantype)
        RefldV  = Canopychar(7,Cantype)
        RefldN  = Canopychar(8,Cantype)
        Cluster = Canopychar(9,Cantype)
!        print*,'cluster',  Cluster
! Extinction coefficients for black leaves for beam (kb) or diffuse (kd)
        Kb = Cluster * 0.5 / Sinbeta 
! (0.5 assumes a spherical leaf angle distribution (0.5 = cos (60 deg))
        Kd = 0.8 * Cluster              
! (0.8 assumes a spherical leaf angle distribution)

        Call CalcExtCoeff(Qbeamv,ScatV,Kb,Kd,ReflbV,KbpV,KdpV,QbAbsV)
        Call CalcExtCoeff(Qbeamn,ScatN,Kb,Kd,ReflbN,KbpN,KdpN,QbAbsn)

        DO i = 1,layers
! LAI depth at this layer
          LAIdepth   = LAI  * Distgauss(i)  
!fraction of leaves that are sunlit
          Sunfrac(i) = EXP(-Kb * LAIdepth)  

          Call CalcRadComponents(Qdiffv , Qbeamv , kdpV,  
     &                          kbpV, kb, scatV, refldV, 
     &                          reflbV, LAIdepth, QdAbsVL, QsAbsVL)

          Call CalcRadComponents(Qdiffn , Qbeamn , kdpN,  
     &                          kbpN, kb, scatN, refldN, 
     &                          reflbN, LAIdepth, QdAbsNL, QsAbsNL)

       ShadePPFD(i) = (QdAbsVL + QsAbsVL) * ConvertShadePPFD/(1 - scatV)
       SunPPFD(i) = ShadePPFD(i) + (QbAbsV* ConvertSunPPFD/(1 - scatV))
          QdAbsV(i) = QdAbsVL
          QsAbsV(i) = QsAbsVL
          QdAbsn(i) = QdAbsNL
          QsAbsn(i) = QsAbsNL
          ShadeQv(i) = QdAbsVL + QsAbsVL
          SunQv(i)   = ShadeQv(i) + QbAbsV 
          ShadeQn(i) = QdAbsNL + QsAbsNL
          SunQn(i)   = ShadeQn(i) + QbAbsn 
        ENDDO

      ELSE                           ! Night time

        QbAbsV  = 0
        QbAbsn   = 0

        Sunfrac(:)   = 0.2
        SunQn(:)     = 0
        ShadeQn(:)   = 0
        SunQv(:)     = 0
        ShadeQv(:)   = 0
        SunPPFD(:)   = 0
        ShadePPFD(:) = 0
        QdAbsV(:)    = 0
        QsAbsV(:)    = 0
        QdAbsn(:)    = 0
        QsAbsn(:)    = 0

      ENDIF
         
      RETURN
      END SUBROUTINE CanopyRad

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   Subroutine CalcExtCoeff
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      SUBROUTINE CalcExtCoeff(Qbeam,scat,kb,kd,reflb,
     &                        kbp,kdp,QbeamAbsorb)

      IMPLICIT NONE
      INTEGER ,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)

      REAL,INTENT(IN) :: Qbeam, scat, Kb, Kd
      REAL,INTENT(OUT) :: Reflb, Kbp, Kdp, QbeamAbsorb

! local variables
      REAL :: P
!-------------------------------------------------------------------
      
      P     = (1 - scat)**0.5
      Reflb = 1 - Exp((-2 * ((1 - P) / (1 + P)) * kb) / (1 + kb))

! Extinction coefficients
      Kbp   = Kb * P
      Kdp   = Kd * P
! Absorbed beam radiation
      QbeamAbsorb = kb * Qbeam * (1 - scat) 

      RETURN
      END SUBROUTINE CalcExtCoeff

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   Subroutine CalcRadComponents
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      SUBROUTINE CalcRadComponents(Qdiff, Qbeam, kdp, kbp, kb,    
     &        scat, refld, reflb, LAIdepth, QdAbs, QsAbs)

      IMPLICIT NONE
      INTEGER ,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)

      REAL,INTENT(IN) :: Qdiff,Qbeam,kdp,kbp,kb,scat,
     &                   refld,reflb,LAIdepth
      REAL,INTENT(OUT) :: QdAbs, QsAbs
!-------------------------------------------------------------------

      QdAbs = Qdiff *   Kdp * (1 - Refld) * Exp(-Kdp * LAIdepth)
      QsAbs = Qbeam * ((Kbp * (1 - Reflb) * Exp(-Kbp * LAIdepth)) - 
     &                  (Kb * (1 - Scat) * Exp(-Kb * LAIdepth)))

      RETURN
      END SUBROUTINE CalcRadComponents

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   Subroutine CanopyEB
!
!   Canopy energy balance model for estimating leaf temperature
!   Code developed by Alex Guenther, based on Goudrian and Laar (1994),
!    Leuning (1997)
!   Initial code 8-99, modified 7-2000 and 12-2001
!
!   Note: i denotes an array containing a vertical profile through the 
!         canopy with 0 
!         (above canopy conditions) plus 1 to number of canopy layers
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      SUBROUTINE CanopyEB(Trate, Layers, Distgauss, Canopychar,    
     &             Cantype, StomataDI, TairK, HumidairPa, Ws,       
     &             SunPPFD, ShadePPFD, SunQv, ShadeQv, SunQn, ShadeQn,  
     &             Sunleaftk, SunleafSH, SunleafLH,                 
     &             SunleafIR, Shadeleaftk, ShadeleafSH,             
     &             ShadeleafLH, ShadeleafIR, NrCha, NrTyp, Ws0,     
     &             TairK0, HumidairPa0)

      IMPLICIT NONE
      INTEGER ,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)

! inputs
      INTEGER,INTENT(IN) :: NrCha, NrTyp, Layers, Cantype
      REAL,INTENT(IN) :: Trate, StomataDI, TairK0, HumidairPa0, Ws0
      REAL,DIMENSION(Layers),INTENT(IN) ::  Distgauss, SunQv,ShadeQv, 
     &            SunQn, ShadeQn, SunPPFD, ShadePPFD
      REAL,DIMENSION(NrCha, NrTyp),INTENT(IN)  :: Canopychar

! outputs
      REAL,DIMENSION(Layers),INTENT(OUT) :: HumidairPa,            
     &       Ws, Sunleaftk, SunleafSH, SunleafLH, SunleafIR, TairK, 
     &       Shadeleaftk, ShadeleafSH, ShadeleafLH, ShadeleafIR

! local variables
      INTEGER :: i
      REAL :: Cdepth, Lwidth, Llength, Cheight, Eps, TranspireType,   
     &         Deltah, UnexposedLeafIRin, ExposedLeafIRin, IRin,IRout
      REAL,DIMENSION(Layers) :: Ldepth, Wsh
!-----------------------------------------------------------------------
         
      Cdepth        = Canopychar(1, Cantype)
      Lwidth        = Canopychar(2, Cantype)
      Llength       = Canopychar(3, Cantype)
      Cheight       = Canopychar(4, Cantype)
      Eps           = Canopychar(10,Cantype)
      TranspireType = Canopychar(11,Cantype)

      IF (TairK0  > 288) THEN
! Pa m-1  (humidity profile for T < 288)
        Deltah =  Canopychar(14, Cantype) / Cheight
      ELSEIF (TairK0  > 278) THEN
        Deltah =(Canopychar(14,Cantype)-((288-TairK0)/10) *             
     &          (Canopychar(14,Cantype)-Canopychar(15,Cantype)))/Cheight
      ELSE
! Pa m-1  (humidity profile for T <278)
        Deltah = Canopychar(15, Cantype) / Cheight
      ENDIF

      Ldepth(:)     = Cdepth * Distgauss(:)
      TairK(:)      = TairK0  + (Trate  * Ldepth(:))      ! check this
      HumidairPa(:) = HumidairPa0  + (Deltah * Ldepth(:))

      Wsh(:) = (Cheight-Ldepth(:)) - (Canopychar(16,Cantype) * Cheight)
      Ws(:)  = (Ws0*LOG(Wsh(:))/LOG(Cheight-Canopychar(16,Cantype)
     &           *Cheight))
      WHERE (Wsh(:) < 0.001) Ws(:) = 0.05

      DO i=1,Layers
        IRin     = UnexposedLeafIRin(TairK(i), Eps)
        ShadeleafIR(i) = 2 * IRin
        SunleafIR(i) = 0.5*ExposedLeafIRin(HumidairPa0,TairK0)+1.5*IRin

      ! Sun
        CALL LeafEB(SunPPFD(i), SunQv(i) + SunQn(i),                    
     &               SunleafIR(i), Eps, TranspireType, Lwidth, Llength, 
     &               TairK(i), HumidairPa(i), Ws(i),                   
     &               Sunleaftk(i), SunleafSH(i),SunleafLH(i),          
     &               IRout,StomataDI )

         SunleafIR(i) = SunleafIR(i) - IRout

      ! Shade
        CALL LeafEB(ShadePPFD(i), ShadeQv(i)+ShadeQn(i),               
     &               ShadeleafIR(i),Eps,TranspireType, Lwidth,Llength,
     &               TairK(i), HumidairPa(i), Ws(i),        
     &               Shadeleaftk(i), ShadeleafSH(i),ShadeleafLH(i),
     &               IRout, StomataDI )

         ShadeleafIR(i) = ShadeleafIR(i) - IRout
      ENDDO

      RETURN
      END SUBROUTINE CanopyEB

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   Subroutine LeafEB
!
!   Leaf energy balance
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      SUBROUTINE LeafEB(PPFD, Q, IRin, Eps, TranspireType,             
     &         Lwidth, Llength, TairK, HumidairPa, Ws, Tleaf,           
     &         SH, LH, IRout, StomataDI)

      IMPLICIT NONE
      INTEGER,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)

      REAL,INTENT(IN) :: Eps, TranspireType, Lwidth, Llength,PPFD, Q,  
     &                    IRin, TairK, HumidairPa, StomataDI, Ws
      REAL,INTENT(OUT) :: IRout, Tleaf, SH, LH

! local variables

      INTEGER :: i
      REAL :: HumidAirKgm3,GHforced,StomRes,IRoutairT,LHV,LatHv,Ws1, 
     &        LHairT,Tdelt,Balance,LeafBLC,LeafH,LeafLE,LeafIRout,   
     &        GH1,SH1,LH1,E1,ConvertHumidityPa2kgm3,ResSC,IRout1,GH
!----------------------------------------------------

      IF (Ws <= 0) THEN
        Ws1 = 0.001
      ELSE
        Ws1 = Ws
      ENDIF

      ! Air vapor density kg m-3
      HumidAirKgm3 = ConvertHumidityPa2kgm3(HumidairPa, TairK)

      ! Heat convection coefficient (W m-2 K-1) for forced convection. 
      ! Nobel page 366
      GHforced = 0.0259 / (0.004 * ((Llength / Ws)**0.5))

      ! Stomatal resistence s m-1
      StomRes  = ResSC(PPFD, stomataDI)

      IRoutairT = LeafIROut(tairK, eps)

      ! Latent heat of vaporization (J Kg-1)
      LatHv = LHV(TairK)

      ! Latent heat flux
      LHairT = LeafLE(TairK,HumidAirKgm3,LatHv,GHforced,StomRes,
     &                TranspireType)

      E1 = (Q + IRin - IRoutairT - LHairT)
      IF (E1 .EQ. 0.) THEN
        E1 = -1.
      ENDIF

      Tdelt = 1.
      Balance = 10.
      DO i = 1, 10
        IF (ABS(Balance) > 2) THEN
          ! Boundary layer conductance
          GH1 = LeafBLC(GHforced, Tdelt, Llength)
          ! Convective heat flux
          SH1 = LeafH(Tdelt, GH1)                
          ! Latent heat of vaporization (J Kg-1)
          LatHv = LHV(TairK + Tdelt)             
          LH = LeafLE(TairK + Tdelt, HumidAirKgm3,      
     &                 LatHv, GH1, StomRes, TranspireType)
          LH1 = LH - LHairT
          IRout  = LeafIROut(TairK + Tdelt, Eps)
          IRout1 = IRout - IRoutairT
          Tdelt  = E1 / ((SH1 + LH1 + IRout1) / Tdelt)
          Balance = Q + IRin - IRout - SH1 - LH
        ENDIF
      ENDDO

      If (Tdelt > 10)  Tdelt = 10
      If (Tdelt < -10) Tdelt = -10

      Tleaf = TairK + Tdelt
      GH    = LeafBLC(GHforced, Tleaf - TairK, Llength)
      SH    = LeafH(Tleaf - TairK, GH)
      LH    = LeafLE(Tleaf, HumidAirKgm3, LatHv,         
     &                GH, StomRes, TranspireType)
      IRout = LeafIROut(Tleaf, Eps)

      RETURN
      END SUBROUTINE LeafEB

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION Calcbeta
!   Calculates the solar zenith angle
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      FUNCTION Calcbeta(Day, Lat, Hour)

      IMPLICIT NONE
      INTEGER,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)

      INTEGER :: Day

      REAL :: Rpi, Hour, Lat, SinDelta,                 
     &          CosDelta, A, B, Sinbeta, Calcbeta
      REAL,PARAMETER :: Pi = 3.14159, Rpi180 = 57.29578
!--------------------------------------------------------------------
      SinDelta = -SIN(0.40907) * COS(6.28 * (Day + 10) / (365))
      CosDelta = (1 - SinDelta**2.)**0.5

      A = SIN(Lat / Rpi180) * SinDelta
      B = COS(Lat / Rpi180) * Cosdelta
      Sinbeta = A + B * COS(2 * Pi * (Hour - 12) / 24)
      Calcbeta = ASIN(Sinbeta) * 57.29578

      END FUNCTION Calcbeta

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION DIstomata
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      FUNCTION DIstomata(DI)

      IMPLICIT NONE
      INTEGER ,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
      REAL :: DI, DIstomata
! > -.5 incipient,  mild or no drought; < -4 extreme drought
      REAL,PARAMETER :: DIhigh = -0.5, DIlow = -5
!--------------------------------------------------------------------

      IF (DI > DIhigh) THEN
        DIstomata = 1  ! no drought
      ELSEIF (DI > DIlow) THEN
! interpolate
        DIstomata = 1 - (0.9 * ((DI - DIhigh) / (DIlow - DIhigh))) 
      ELSE
        DIstomata = 0  ! Maximum drought, maximum stomatal resistance
      ENDIF

      END FUNCTION DIstomata

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION CalcEccentricity
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      FUNCTION CalcEccentricity(Day)

      IMPLICIT NONE
      INTEGER :: Day
      INTEGER,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
      REAL :: CalcEccentricity
!--------------------------------------------------------------------

      CalcEccentricity = 1 + 0.033 * COS(2*3.14*(Day-10)/365)
    
      END FUNCTION CalcEccentricity

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION UnexposedLeafIRin
!
!   Calculate IR into leaf that is not exposed to the sky
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

       FUNCTION UnexposedLeafIRin(Tk, Eps)

       IMPLICIT NONE
       INTEGER ,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
       REAL :: Eps, Tk, UnexposedLeafIRin
! Stefan-boltzman constant  W m-2 K-4
       REAL,PARAMETER :: Sb = 0.0000000567 
!--------------------------------------------------------------------

       UnexposedLeafIRin = Eps * Sb * (Tk**4.) 

       END FUNCTION UnexposedLeafIRin

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION ExposedLeafIRin
!
!   Calculate IR into leaf that is exposed to the sky
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

       FUNCTION ExposedLeafIRin(HumidPa, Tk)

       IMPLICIT NONE
       INTEGER ,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
       REAL :: Tk, HumidPa, EmissAtm, ExposedLeafIRin
! Stefan-boltzman constant  W m-2 K-4
       REAL,PARAMETER :: Sb = 0.0000000567 
!--------------------------------------------------------------------

! Apparent atmospheric emissivity for clear skies: 
! function of water vapor pressure (Pa) 
! and ambient Temperature (K) based on Brutsaert(1975) 
! referenced in Leuning (1997)

        EmissAtm        = 0.642 * (HumidPa / Tk)**(1./7.)
        ExposedLeafIRin = EmissAtm * Sb * (Tk**4.)

        END FUNCTION ExposedLeafIRin

!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION WaterVapPres
!
!   Convert water mixing ratio (kg/kg) to water vapor pressure 
!   (Pa or Kpa depending on units of input )
!   Mixing ratio (kg/kg), temp (C), pressure (KPa)
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      FUNCTION WaterVapPres(Dens, Pres, WaterAirRatio)

      IMPLICIT NONE
      INTEGER ,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
      REAL :: Dens, Pres, WaterVapPres, WaterAirRatio
!--------------------------------------------------------------------

      WaterVapPres = (Dens / (Dens + WaterAirRatio)) * Pres

      END FUNCTION WaterVapPres

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION Stability
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      FUNCTION Stability(Canopychar, Cantype, Solar, NrCha, NrTyp)

      IMPLICIT NONE
      INTEGER :: Cantype, NrCha, NrTyp
      INTEGER ,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
      REAL :: Solar, Trateboundary, Stability
      REAL,DIMENSION(NrCha, NrTyp)  :: Canopychar
!--------------------------------------------------------------------

      Trateboundary = 500

      IF (Solar > Trateboundary) THEN
            ! Daytime temperature lapse rate
        Stability = Canopychar(12, Cantype)
      ELSEIF (Solar > 0) THEN
        Stability = Canopychar(12, Cantype) -                      
     &             ((Trateboundary - Solar) / Trateboundary) *     
     &    (Canopychar(12, Cantype) - Canopychar(13, Cantype))
       ELSE
            ! Nightime temperature lapse rate
         Stability = Canopychar(13, Cantype)
       ENDIF

       END FUNCTION Stability

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION ConvertHumidityPa2kgm3
!
!   Saturation vapor density  (kg/m3)
!
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      FUNCTION ConvertHumidityPa2kgm3(Pa, Tk)

      IMPLICIT NONE
      INTEGER,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
      REAL :: ConvertHumidityPa2kgm3, Pa, Tk
!--------------------------------------------------------------------

      ConvertHumidityPa2kgm3 = 0.002165 * Pa / Tk

      END FUNCTION ConvertHumidityPa2kgm3

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION ResSC
!
!   Leaf stomatal cond. resistance s m-1
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      FUNCTION ResSC(Par, StomataDI)

      IMPLICIT NONE
      INTEGER,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
      REAL :: Par, StomataDI, SCadj, ResSC
!--------------------------------------------------------------------

      SCadj = StomataDI * ((0.0027 * 1.066 * Par) /    
     &        ((1 + 0.0027 * 0.0027 * Par**2.)**0.5))

      IF (SCadj < 0.1) THEN
        ResSC = 2000
      ELSE
        ResSC = 200 / SCadj
      ENDIF

      END FUNCTION ResSC

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION LeafIROut
!
!   IR thermal radiation energy output by leaf
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      FUNCTION LeafIROut(Tleaf, Eps)

      IMPLICIT NONE
      INTEGER ,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
      REAL :: Tleaf, Eps, LeafIROut
! Stefan-boltzman constant  W m-2 K-4
      REAL,PARAMETER :: Sb = 0.0000000567 
!--------------------------------------------------------------------

      LeafIROut = Eps * Sb * (2 * (Tleaf**4.))

      END FUNCTION LeafIROut

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION LHV
!
!   Latent Heat of vaporization(J Kg-1) from Stull p641
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      FUNCTION LHV(Tk)

      IMPLICIT NONE
      INTEGER ,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
      REAL :: Tk, LHV
!--------------------------------------------------------------------

      LHV = 2501000 - (2370 * (Tk - 273))

      END FUNCTION LHV

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION LeafLE
!
!   Latent energy term in Energy balance
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      FUNCTION LeafLE(Tleaf, Ambvap, LatHv, GH, StomRes, TranspireType)

      IMPLICIT NONE
      INTEGER ,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
      REAL :: Tleaf, Ambvap, LatHv, GH, StomRes,                      
     &         TranspireType, SvdTk,LeafRes, Vapdeficit, LeafLE, LE
!--------------------------------------------------------------------

      LeafRes    = (1 / (1.075 * (GH / 1231))) + StomRes
      Vapdeficit = (SvdTk(Tleaf) - Ambvap)

! Latent heat of vap (J Kg-1) * vap deficit(Kg m-3) / 
!                 leaf resistence (s m-1)
      LE = TranspireType * (1 / LeafRes) * LatHv * Vapdeficit
      IF (LE < 0) THEN
        LeafLE = 0
      ELSE
        LeafLE = LE
      ENDIF

      END FUNCTION  LeafLE

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION LeafBLC
!
!   Boundary layer conductance
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      FUNCTION LeafBLC(GHforced, Tdelta, Llength)

      IMPLICIT NONE
      INTEGER ,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
      REAL :: GHforced, Tdelta, Llength, Ghfree, LeafBLC
!--------------------------------------------------------------------

! This is based on Leuning 1995 p.1198 except using molecular 
! conductivity (.00253 W m-1 K-1 Stull p 640) instead of molecular
! diffusivity so that you end up with a heat convection coefficient 
! (W m-2 K-1) instead of a conductance for free convection

      IF (Tdelta >= 0) THEN
         GhFree = 0.5 * 0.00253 * ((160000000 * Tdelta /            
     &             (Llength**3.))**0.25) / Llength
      ELSE
        GhFree = 0
      ENDIF
      LeafBLC = GHforced + GhFree

      END FUNCTION LeafBLC

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION LeafH
!
!   Convective energy term in Energy balance (W m-2 heat flux from 
!      both sides of leaf)
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      FUNCTION LeafH(Tdelta, GH)
  
      IMPLICIT NONE
      INTEGER,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
      REAL :: Tdelta, GH, LeafH
!--------------------------------------------------------------------

! 2 sides X conductance X Temperature gradient
      LeafH = 2 * GH * Tdelta
  
      END FUNCTION LeafH

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION SvdTk
!
!   Saturation vapor density  (kg/m3)
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      Function SvdTk(Tk)

      IMPLICIT NONE
      INTEGER ,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
      REAL :: Tk, Svp, SvdTk
!--------------------------------------------------------------------

! Saturation vapor pressure (millibars)
      Svp = 10**((-2937.4 / Tk) - (4.9283 * LOG10(Tk)) + 23.5518)  
      SvdTk = 0.2165 * Svp / Tk

      END FUNCTION  SvdTk

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION Ea1t99
!
!   Temperature dependence activity factor for emission type 1 
!          (e.g. isoprene, MBO)
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      FUNCTION Ea1t99(T1, T24, T240,SPC_NAME)

      IMPLICIT NONE
      INCLUDE 'CONST_CANOPY.EXT'
      INTEGER,EXTERNAL :: INDEX1
        REAL,PARAMETER :: Ctm2 = 230
      REAL :: T1, T24, T240, Ea1t99, Topt, X, Eopt
!--------------------------------------------------------------------
       CHARACTER*16,INTENT(IN) :: SPC_NAME
       INTEGER :: SPCNUM
         SPCNUM = INDEX1(SPC_NAME,N_MGN_SPC,MGN_SPC)
      IF (T1 < 260) THEN
        Ea1t99 = 0
      ELSE
        ! Energy of activation and deactivation

        ! Temperature at which maximum emission occurs
        Topt = 312.5 + 0.6 * (T240 - 297)
        X    = ((1 / Topt) - (1 / T1)) / 0.00831

        ! Maximum emission (relative to emission at 30 C)
        Eopt = CLeo(SPCNUM) * EXP(0.05 * (T24 - 297)) *
     &         Exp(0.05*(T240-297))          
        Ea1t99 = Eopt * Ctm2 * Exp(Ctm1(SPCNUM) * X) /           
     &       (Ctm2 - Ctm1(SPCNUM) * (1 - EXP(Ctm2 * X)))
!       print*,'SPCNUM---',SPCNUM,CLeo(SPCNUM),Ctm1(SPCNUM)
      ENDIF
      

      END FUNCTION  Ea1t99

!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION Ea1pp
!
! pstd = 200 for sun leaves and 50 for shade leaves 
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      FUNCTION Ea1p99(PPFD1, PPFD24, PPFD240, PSTD)

      IMPLICIT NONE
      REAL :: PPFD1, PPFD24, PPFD240, PSTD, Alpha, C1, Ea1p99
!--------------------------------------------------------------------

      IF (PPFD240 < 0.01) THEN
        Ea1p99 = 0
      ELSE
        Alpha  = 0.004 - 0.0005 * LOG(PPFD240)
        C1     = 0.0468 * EXP(0.0005 * (PPFD24 - PSTD))   
     &          * (PPFD240 ** 0.6)
        Ea1p99 = (Alpha * C1 * PPFD1) /                   
     &          ((1 + Alpha**2. * PPFD1**2.)**0.5)
      ENDIF 


      END FUNCTION  Ea1p99

!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION Ealti99
!
! calculate light indepent algorithms
! coded by Xuemei Wang 05 Nov. 2007
!--   GAMMA_TLI =  exp[BETA*(T-Ts)]
!           where BETA   = temperature dependent parameter
!                 Ts     = standard temperature (normally 303K, 30C)
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      FUNCTION Ealti99(spcnam, temp)

      IMPLICIT NONE
      INCLUDE 'EACO.EXT'

      INTEGER :: spcnum                             ! Species number
      INTEGER,EXTERNAL :: index1
      CHARACTER(LEN=16) :: spcnam
      REAL :: temp, Ealti99
!     REAL,PARAMETER :: Ts = 303.0
      REAL,PARAMETER :: Ts = 303.15 
!--------------------------------------------------------------------

      spcnum = index1(spcnam,n_tdf_spc,tdf_spc)
      Ealti99 = exp( tdf_prm(spcnum)*(temp-Ts) )


      END FUNCTION Ealti99
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
