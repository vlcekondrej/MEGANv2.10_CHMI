      PROGRAM EMPROC

!***********************************************************************
!   This program computes biogenic emission using input emission 
!   capacity maps and MCIP output variables.
!   The emission capacity map (inpname) are gridded in netCDF-IOAPI format
!   with all the daily average ppfd and daily average temperature.
!
!   Note: The projection and input grids of the two files must be 
!   identical.
!
!   Requirement:
!      Requires libnetcdf.a and libioapi.a to compile
!
!      setenv EFMAPS  <input netCDF-ioapi emission capacity map>
!      setenv PFTS16  <input netCDF-ioapi PFT data>
!      setenv LAIS46  <input netCDF-ioapi LAI data>
!      setenv MGNMET  <input netCDF-ioapi MEGAN met>
!      setenv EROUT <megan emission output>
!      setenv SDATE <start date>
!      setenv STIME <start time>
!      setenv RLENG <run length>
!
!      If SDATE,STIME and RLENG are not set, use values from MET input
!
!   CALL:
!      CHECKMEM
!      CKGDIOAPI2D
!      MODULE GAMMA_ETC
!         GAMMA_LAI
!         GAMMA_P
!         GAMMA_TLD
!         GAMMA_TLI
!         GAMMA_A
!         GAMMA_S
!
!   History:
!   Created by Jack Chen 11/04
!   Modified by Tan 11/21/06 for MEGAN v2.0
!   Modified by Xuemei Wang 11/04/2007 for MEGAN2.1
!   Modified by Julia Lee-Taylor 03/18/2008 for MEGAN2.1
!   Modified by Xuemei Wang 09/30/2008 for MEGAN2.1
!   Modified by Tan         07/28/2011 for MEGAN2.1
!
!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Scientific algorithm
!
!             Emission = [EF][GAMMA][RHO]
!           where [EF]    = emission factor (ug/m2h)
!                 [GAMMA] = emission activity factor (non-dimension)
!                 [RHO]   = production and loss within plant canopies
!                           (non-dimensional)
!                 Assumption: [RHO] = 1  (11/27/06) (See PDT_LOT_CP.EXT)
!
!             GAMMA  = [GAMMA_CE][GAMMA_age][GAMMA_SM]
!           where [GAMMA_CE]  = canopy correction factor
!                 [GAMMA_AGE] = leaf age correction factor
!                 [GAMMA_SM]  = soil moisture correction factor
!                 Assumption: [GAMMA_SM]  = 1  (11/27/06)

!             GAMMA_CE = [GAMMA_LAI][GAMMA_P][GAMMA_T]
!           where [GAMMA_LAI] = leaf area index factor
!                 [GAMMA_P]   = ppfd emission activity factor
!                 [GAMMA_T]   = temperature response factor
!
!             Emission = [EF][GAMMA_LAI][GAMMA_P][GAMMA_T][GAMMA_AGE]
!        Derivation:
!             Emission = [EF][GAMMA](1-LDF) + [EF][GAMMA][LDF][GAMMA_P]
!             Emission = [EF][GAMMA]{ (1-LDF) + [LDF][GAMMA_P] }
!             Emission = [EF][GAMMA]{ (1-LDF) + [LDF][GAMMA_P] }
!           where LDF = light dependent function (non-dimension)
!                               (See LD_FCT.EXT)
!
!        Final Equation
!             Emission = [EF][GAMMA_LAI][GAMMA_AGE]*
!    { (1-LDF)[GAMMA_TLI] + [LDF][GAMMA_P][GAMMA_TLD] }  !for MEGAN2.1
!         where GAMMA_TLI is light independent
!               GAMMA_TLD is light dependent  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      USE GAMMA_ETC                 ! Module containing gamma functions

      IMPLICIT NONE

!...  INCLUDES:
      INCLUDE 'PARMS3.EXT'          !  I/O API parameters
      INCLUDE 'IODECL3.EXT'         !  I/O API function declarations
      INCLUDE 'FDESC3.EXT'          !  I/O API file description data structures
      INCLUDE 'EACO.EXT'            ! emissions factor type definitions
      INCLUDE 'CONST_MEGAN.EXT'

!...  EXTERNAL FUNCTIONS and their descriptions:
      INTEGER, EXTERNAL   ::  FINDC, ENVINT
      LOGICAL, EXTERNAL   ::  CKGDIOAPI2D, ENVYN
      INTEGER      PROMPTFFILE, STR2INT
      REAL         STR2REAL
      INTEGER      INDEX1
      EXTERNAL     PROMPTFFILE, STR2INT
      EXTERNAL     STR2REAL
      EXTERNAL     INDEX1
      LOGICAL      DSCGRID
      EXTERNAL     DSCGRID

!...  Program I/O files: From run script
! Program name
      CHARACTER*16  :: PROGNAME = 'EMPROC'
! Netcdf file
      CHARACTER*16  :: EFMAPS = 'EFMAPS'     ! Output EF file logical name
      CHARACTER*16  :: PFTS16 = 'PFTS16'     ! Output PFT file logical name
      CHARACTER*16  :: LAIS46 = 'LAIS46'     ! Output LAI file logical name
! Met files
      CHARACTER*16  :: MGNMET = 'MGNMET'     ! Met file logical name
! Output file
      CHARACTER*16  :: MGNERS = 'MGNERS'     ! Output file logical name

!...  Parameters for file units
      INTEGER  LOGDEV                      ! Logfile unit number

!...  External parameters
! From run script
      INTEGER       SDATE          ! Start date YYYYDDD
      INTEGER       STIME          ! Start time HHMMSS
      INTEGER       RLENG          ! Run length HHMMSS

      LOGICAL       ONLN_DT        ! Daily average, temperature
      LOGICAL       ONLN_DS        ! Daily average, solar radiation 

! I/O API file parameters
      INTEGER :: JDATE        ! Date YYYYDDD from inpname
      INTEGER :: JTIME        ! Time HHMMSS from inpname
      INTEGER :: NCOLS        ! Number of columns
      INTEGER :: NROWS        ! Number of rows
      INTEGER :: NLAYS        ! Number of vertical layers
      INTEGER :: MXREC        ! Total number of timesteps
      INTEGER :: TSTEP        ! Time step

!...  Internal parameters
! Internal parameters (status and buffer)
      INTEGER       IOS                    ! i/o status
      CHARACTER*256 MESG                   ! message buffer

! Parameter for output species
      INTEGER, PARAMETER :: NEMIS = N_MGN_SPC 
                          ! number of output emission variables
                          ! number of MEGAN species

! Local variables and their descriptions:
      CHARACTER*16  :: GDNAM
      CHARACTER*16  :: CNAME        ! Coord name

      REAL          GAREA           ! Area in one grid (metre^2)
      INTEGER       AVEBY           ! Divider for daily average

      REAL          LDF             ! Light dependent factor
      REAL          RHO             ! Production and loss within canopy

      REAL, ALLOCATABLE :: ER( :,: )      ! Output emission buffer
      REAL, ALLOCATABLE :: NON_DIMGARMA (:,:,:)
      REAL, ALLOCATABLE :: LAT( :,: )     ! Latitude of grid cell
      REAL, ALLOCATABLE :: LONG( :,: )    ! Longitude of grid cell

      REAL, ALLOCATABLE :: LAIp( :,: )    ! Previous monthly LAI
      REAL, ALLOCATABLE :: LAIc( :,: )    ! Current monthly LAI

      REAL, ALLOCATABLE :: TEMP( :,: )    ! Temperature (K)
      REAL, ALLOCATABLE :: PPFD( :,: )    ! Calculated PAR (umol/m2.s)

      REAL, ALLOCATABLE :: D_PPFD( :,: )  ! Daily PAR (umol/m2.s)
      REAL, ALLOCATABLE :: D_TEMP( :,: )  ! Daily temperature (K)
      REAL, ALLOCATABLE :: GAM_LHT( :,: ) ! LAI correction factor
      REAL, ALLOCATABLE :: GAM_AGE( :,: ) ! leaf age correction factor
      REAL, ALLOCATABLE :: GAM_SMT( :,: ) ! Soil moisture correction factor
      REAL, ALLOCATABLE :: WIND( :,: )
      REAL, ALLOCATABLE :: PRES( :,: )
      REAL, ALLOCATABLE :: DI( :,: )
      REAL, ALLOCATABLE :: QV( :,: )
      REAL, ALLOCATABLE :: PRECADJ( :,: )

      INTEGER :: I_PFT
      INTEGER :: LAIp_DY, LAIp_HR, LAIc_DY, LAIc_HR
      INTEGER :: MXPFT, MXLAI

! Number of LAT, LONG, and PFT factor variables
      REAL, ALLOCATABLE :: PFTF( :, :, : )   ! PFT factor array
      REAL, ALLOCATABLE :: GAM_TLD(:,:)
      REAL, ALLOCATABLE :: GAM_TLI(:,:)
      REAL :: ADJUST_FACTOR_LD,ADJUST_FACTOR_LI
      REAL :: GAMMA_TD,GAMMA_TI, TotalPFT
        
      ! loop indices
      INTEGER ::  T, S, I, J , K, N     ! Counters
      INTEGER ::  NMAP            ! Index
      INTEGER ::  IDATE           ! Looping
      INTEGER ::  ITIME           ! Looping
! SOIL NOx
      LOGICAL :: LSOIL = .TRUE.
      INTEGER, ALLOCATABLE :: SLTYP( :,: )
      REAL, ALLOCATABLE    :: RSTYP( :,: )
      REAL, ALLOCATABLE    :: SOILM( :,: )
      REAL, ALLOCATABLE    :: SOILT( :,: )
      REAL, ALLOCATABLE    :: CFNO( :,: )
      REAL, ALLOCATABLE    :: CFNOG( :,: )

!...  Constants
! Length of the time step  (days)
      INTEGER, PARAMETER :: TSTLEN = 30

! parameter for unit conversion
      REAL,PARAMETER :: ug2tonne = 1E-12  ! convert microgram to metric tonne
      REAL,PARAMETER :: hr2sec = 3600     ! convert hr to second
      REAL,PARAMETER :: ug2g = 1E-6       ! convert microgram to gram


!***********************************************************************

!--=====================================================================
!...  Begin program
!--=====================================================================

!-----------------------------------------------------------------------
!.....1) File set up and assign I/O parameters
!-----------------------------------------------------------------------
!...  Initialize log file unit
      LOGDEV = INIT3()
           !  Now I/O API is set up, and LOGUNIT is the unit number
           !  for the log file (or it 6 for st'd output).

!...  Check logical variables
      MESG = 'Use daily average temperature'
      ONLN_DT    = ENVYN ( 'ONLN_DT', MESG, .TRUE., IOS )

      MESG = 'Use daily average PAR'
      ONLN_DS    = ENVYN ( 'ONLN_DS', MESG, .TRUE., IOS )

!...  Get input parameters from run script
      CALL ENVSTR( 'GDNAM3D', MESG, 'ASACA36km', GDNAM, IOS )
      IF( .NOT. DSCGRID( GDNAM, CNAME, GDTYP3D,
     &              P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,
     &              XORIG3D, YORIG3D, XCELL3D, YCELL3D,
     &              NCOLS3D, NROWS3D, NTHIK3D ) ) THEN
         MESG = 'Could not get grid description.'
         CALL M3EXIT ( PROGNAME, 0, 0, MESG, 2 )
      ENDIF

!...  Open files
      WRITE(MESG,1030) 'Checking up files',0,0,0
      CALL M3MESG( MESG )
      IF ( .NOT. OPEN3( PFTS16, FSREAD3, PROGNAME ) ) THEN
         CALL NAMEVAL (PFTS16, MESG)  ! get input file name and path
         MESG = 'Could not open file '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      ! Check grid
      IF ( .NOT. FILCHK3 ( PFTS16,
     &              GRDDED3, NCOLS3D, NROWS3D, 1, NTHIK3D))  THEN
         MESG = 'PFTS16 has differenet grid definition'
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      IF ( .NOT. DESC3( PFTS16 ) ) THEN
         CALL NAMEVAL (PFTS16, MESG)  ! get input file name and path
         MESG = 'Could not get description of '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      MXPFT = MXREC3D


      IF ( .NOT. OPEN3( LAIS46, FSREAD3, PROGNAME ) ) THEN
         CALL NAMEVAL (LAIS46, MESG)  ! get input file name and path
         MESG = 'Could not open file '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      ! Check grid
      IF ( .NOT. FILCHK3 ( LAIS46,
     &              GRDDED3, NCOLS3D, NROWS3D, 1, NTHIK3D))  THEN
         MESG = 'LAIS46 has differenet grid definition'
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      IF ( .NOT. DESC3( LAIS46 ) ) THEN
         CALL NAMEVAL (LAIS46, MESG)  ! get input file name and path
         MESG = 'Could not get description of '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      MXLAI = MXREC3D

      IF ( .NOT. OPEN3( EFMAPS, FSREAD3, PROGNAME ) ) THEN
         CALL NAMEVAL (EFMAPS, MESG)  ! get input file name and path
         MESG = 'Could not open file '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      ! Check grid
      IF ( .NOT. FILCHK3 ( EFMAPS,
     &              GRDDED3, NCOLS3D, NROWS3D, 1, NTHIK3D))  THEN
         MESG = 'EFMAPS has differenet grid definition'
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      IF ( .NOT. DESC3( EFMAPS ) ) THEN
         CALL NAMEVAL (EFMAPS, MESG)  ! get input file name and path
         MESG = 'Could not get description of '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF

      IF ( .NOT. OPEN3( MGNMET, FSREAD3, PROGNAME ) ) THEN
         CALL NAMEVAL (MGNMET, MESG)  ! get input file name and path
         MESG = 'Could not open file '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      ! Check grid
      IF ( .NOT. FILCHK3 ( MGNMET,
     &              GRDDED3, NCOLS3D, NROWS3D, 1, NTHIK3D))  THEN
         MESG = 'MGNMET has differenet grid definition'
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      IF ( .NOT. DESC3( MGNMET ) ) THEN
         CALL NAMEVAL (MGNMET, MESG)  ! get input file name and path
         MESG = 'Could not get description of '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF


      NCOLS = NCOLS3D
      NROWS = NROWS3D
      TSTEP = TSTEP3D
      GAREA = XCELL3D * YCELL3D

!...  Get input parameters from run script
      MESG = 'Model start date (YYYYDDD)'
      SDATE = ENVINT( 'SDATE', MESG, JDATE, IOS )

      MESG = 'Model start time (HHMMSS)'
      STIME = ENVINT( 'STIME', MESG, JTIME, IOS )

      MESG = 'Model run length (HHMMSS)'
      RLENG = ENVINT( 'RLENG', MESG, MXREC*10000, IOS )

!...  Check start date, start time, end date, end time in MGNMET
      WRITE(MESG,1030) 'Checking up MGNMET',0,0,0
      CALL M3MESG( MESG )
      IDATE = SDATE; ITIME = STIME
      IF ( .NOT. CHECK3( MGNMET, 'TEMP2', IDATE, ITIME ) ) THEN  
         MESG = 'Starting time not on met file'
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      CALL NEXTIME ( IDATE, ITIME, RLENG-10000 )
      IF ( .NOT. CHECK3( MGNMET, 'TEMP2', IDATE, ITIME ) ) THEN
         MESG = 'Ending time not on met file'
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      IDATE = SDATE; ITIME = STIME

!...  Set output parameters that are different from met file and open file 
      SDATE3D = SDATE                ! From run-script
      STIME3D = STIME                ! From run-script
      MXREC3D = RLENG / 10000        ! From run-script
      MXREC   = MXREC3D
      NVARS3D = NEMIS + 7

      DO S = 1, NEMIS
         VNAME3D(S) = TRIM( MGN_SPC( S ) )
         VDESC3D(s) = 'Environmental activity factor for '//
     &                TRIM( MGN_SPC(S) )
         UNITS3D(s) = 'Non-Dimension '
         VTYPE3D(s) = M3REAL
!         print*,'VNAME=', vname3d(s),VDESC3d(s),UNITS3d(s)
      ENDDO

      VNAME3D(NEMIS+1) = 'D_TEMP'
      UNITS3D(NEMIS+1) = 'K'
      VTYPE3D(NEMIS+1) = M3REAL
      VDESC3D(NEMIS+1) = 'Variable  '//'K'

      VNAME3D(NEMIS+2) = 'D_PPFD'
      UNITS3D(NEMIS+2) = 'umol/m2.s'
      VTYPE3D(NEMIS+2) = M3REAL
      VDESC3D(NEMIS+2) = 'Variable  '//'umol/m2.s'

      VNAME3D(NEMIS+3) = 'LAT'
      UNITS3D(NEMIS+3) = ' '
      VTYPE3D(NEMIS+3) = M3REAL
      VDESC3D(NEMIS+3) = ' '

      VNAME3D(NEMIS+4) = 'LONG'
      UNITS3D(NEMIS+4) = ' '
      VTYPE3D(NEMIS+4) = M3REAL
      VDESC3D(NEMIS+4) = ' '

      VNAME3D(NEMIS+5) = 'CFNO'
      UNITS3D(NEMIS+5) = ' '
      VTYPE3D(NEMIS+5) = M3REAL
      VDESC3D(NEMIS+5) = ' '

      VNAME3D(NEMIS+6) = 'CFNOG'
      UNITS3D(NEMIS+6) = ' '
      VTYPE3D(NEMIS+6) = M3REAL
      VDESC3D(NEMIS+6) = ' '

      VNAME3D(NEMIS+7) = 'SLTYP'
      UNITS3D(NEMIS+7) = ' '
      VTYPE3D(NEMIS+7) = M3INT
      VDESC3D(NEMIS+7) = ' '

      CALL NAMEVAL (EFMAPS, MESG)  ! get input file name and path
      FDESC3D(3) = 'Input EFMAPS file: '//TRIM(MESG)

      CALL NAMEVAL (MGNMET, MESG)  ! get input file name and path
      FDESC3D(4) = 'Input MGNMET file: '//TRIM(MESG)
      

      IF ( .NOT. OPEN3( MGNERS, FSCREA3, PROGNAME ) ) THEN
         CALL NAMEVAL (MGNMET, MESG)  ! get input file name and path
         MESG = 'Could not open file '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
!-----------------------------------------------------------------------
!.....2) Process emission rates
!-----------------------------------------------------------------------
!...  Allocate memory
      ALLOCATE ( ER    ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'ER',     PROGNAME )
      ALLOCATE ( LAT   ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'LAT',    PROGNAME )
      ALLOCATE ( LONG  ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'LONG',   PROGNAME )
      ALLOCATE ( LAIp  ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'LAIp',   PROGNAME )
      ALLOCATE ( LAIc  ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'LAIc',   PROGNAME )
      ALLOCATE ( D_PPFD( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'D_PPFD', PROGNAME )
      ALLOCATE ( D_TEMP( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'D_TEMP', PROGNAME )
      ALLOCATE ( PPFD  ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'PPFD',   PROGNAME )
      ALLOCATE ( TEMP  ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'TEMP',   PROGNAME )
      ALLOCATE ( WIND  ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'WIND',   PROGNAME )
      ALLOCATE ( PRES  ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'PRES',   PROGNAME )
      ALLOCATE ( DI    ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'DI',     PROGNAME )
      ALLOCATE ( SOILT ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'SOILT',  PROGNAME )
      ALLOCATE ( QV    ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'QV',     PROGNAME )
      ALLOCATE ( SLTYP ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'SLTYP',  PROGNAME )
      ALLOCATE ( RSTYP ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'RSTYP',  PROGNAME )
      ALLOCATE ( SOILM ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'SOILM',  PROGNAME )
      ALLOCATE ( CFNO  ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'CFNO',   PROGNAME )
      ALLOCATE ( CFNOG ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'CFNOG',  PROGNAME )

      ALLOCATE ( GAM_LHT( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM ( IOS, 'GAM_LHT', PROGNAME )
      ALLOCATE ( GAM_AGE( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM ( IOS, 'GAM_AGE', PROGNAME )
      ALLOCATE ( GAM_SMT( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM ( IOS, 'GAM_SMT', PROGNAME )
      ALLOCATE ( GAM_TLD ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM ( IOS, 'GAM_TLD', PROGNAME )
      ALLOCATE ( GAM_TLI ( NCOLS, NROWS ), STAT = IOS )        
      CALL CHECKMEM ( IOS, 'GAM_TLI', PROGNAME )

      ALLOCATE ( PRECADJ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM     ( IOS, 'PRECADJ', PROGNAME )
      ALLOCATE ( PFTF( NRTYP, NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM ( IOS, 'PFTF', PROGNAME )

      ALLOCATE ( NON_DIMGARMA (NEMIS, NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM ( IOS, 'NON_DIMGARMA', PROGNAME )

!...  Read EFMAPS
   
      IF ( .NOT. READ3(EFMAPS,'LAT',1,1,1,LAT)) THEN
         MESG = 'Error reading LAT'
         CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
      ENDIF

      IF ( .NOT. READ3(EFMAPS,'LONG',1,1,1,LONG)) THEN
         MESG = 'Error reading LONG'
         CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
      ENDIF

!...  Read PFTS16
      DO N = 1, MXPFT
         IF ( .NOT. READ3(PFTS16,'PFTS',1,0,(N-1)*10000,PFTF(N,:,:))) THEN
            MESG = 'Error reading PFTS'
            CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
         ENDIF
      ENDDO

!...  Choosing temperature and PAR
      IF ( ONLN_DT ) THEN
      ! Use  daily temperature
        AVEBY = MIN(24,MXREC)
        ! Start the loop over the time period
        IDATE = SDATE
        ITIME = STIME
        D_TEMP = 0.0
        DO T = 1, AVEBY
          IF ( .NOT. READ3(MGNMET,'TEMP2',1,IDATE,ITIME,TEMP)) THEN
            MESG = 'Error reading TEMP2'
            CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
          ENDIF
          D_TEMP = D_TEMP + TEMP
          CALL NEXTIME( IDATE, ITIME, TSTEP )
        ENDDO
        D_TEMP = D_TEMP/AVEBY
      ELSE
        IF ( .NOT. READ3(EFMAPS,'D_TEMP',1,1,1,D_TEMP)) THEN
          MESG = 'Error reading D_TEMP'
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF
      ENDIF

      IF ( ONLN_DS ) THEN
      ! Use daily PAR
        AVEBY = MIN(24,MXREC)
        ! Start the loop over the time period
        IDATE = SDATE
        ITIME = STIME
        D_PPFD = 0.0
        DO T = 1, AVEBY
          IF ( .NOT. READ3(MGNMET,'PAR',1,IDATE,ITIME,PPFD)) THEN
            MESG = 'Error reading PAR'
            CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
          ENDIF
          !PPFD = PPFD * 4.766
          PPFD = PPFD * 4.5
          D_PPFD = D_PPFD + PPFD
          CALL NEXTIME( IDATE, ITIME, TSTEP )
        ENDDO
        D_PPFD = D_PPFD/AVEBY
      ELSE
        IF ( .NOT. READ3(EFMAPS,'D_SRAD',1,1,1,D_PPFD)) THEN
          MESG = 'Error reading D_SRAD'
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF
!       ppfd: srad - short wave from sun (W/m2)
!       assuming 4.766 (umol m-2 s-1) per (W m-2)
!       assume 1/2 of srad is in 400-700nm band
        !D_PPFD = D_PPFD * 4.766 * 0.5
        D_PPFD = D_PPFD * 4.5 * 0.5
      ENDIF


!...  Start the loop over the time period
      IDATE = SDATE
      ITIME = STIME
      DO T = 1, MXREC
        WRITE(MESG,1030) 'Processing: ',T,IDATE,ITIME
        CALL M3MESG( MESG )
!...  Initialize hourly variables
        TEMP = 0.
        PPFD = 0.
        LAIp = 0.
        LAIc = 0.
         
        IF ( .NOT. READ3(MGNMET,'TEMP2',  ALLAYS3,IDATE,ITIME,TEMP)) THEN
          MESG = 'Error reading temperature'
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

        IF ( .NOT. READ3(MGNMET,'PAR',   ALLAYS3,IDATE,ITIME,PPFD)) THEN
          MESG = 'Error reading PAR'
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF
        !PPFD = PPFD * 4.766
        PPFD = PPFD * 4.5
 
        IF( .NOT. READ3(MGNMET,'WINDSPD',ALLAYS3,IDATE,ITIME,WIND)) THEN
          MESG = 'Error reading wind speed'
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

        IF ( .NOT. READ3(MGNMET,'PRES',  ALLAYS3,IDATE,ITIME,PRES)) THEN
          MESG = 'Error reading pressure'
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF  
 
        IF ( .NOT. READ3(MGNMET,'QV',    ALLAYS3,IDATE,ITIME,QV)) THEN
          MESG = 'Error reading QV'
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

        IF ( .NOT. READ3(MGNMET,'PREC_ADJ',ALLAYS3,IDATE,ITIME,PRECADJ)) THEN
          MESG = 'Error reading precipitation adjustment'
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

        ! Find LAIp and LAIc from date     
        CALL FINDLAI(IDATE,MXLAI,LAIp_DY,LAIp_HR,LAIc_DY,LAIc_HR)
        WRITE(MESG,1020) 'Found LAI current period for YYYYJJJ : ',
     &                    IDATE,LAIc_DY,LAIc_HR
        CALL M3MESG( MESG )
        WRITE(MESG,1020) 'Found LAI previous period for YYYYJJJ : ',
     &                    IDATE,LAIp_DY,LAIp_HR
        CALL M3MESG( MESG )
        IF ( .NOT. READ3(LAIS46,'LAIS',ALLAYS3,LAIp_DY,LAIp_HR,LAIp)) THEN
          MESG = 'Error reading precipitation adjustment'
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF
        IF ( .NOT. READ3(LAIS46,'LAIS',ALLAYS3,LAIc_DY,LAIc_HR,LAIc)) THEN
          MESG = 'Error reading precipitation adjustment'
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

        ! Go over all the chemical species
        DO S = 1, NEMIS
          ! Initialize variables
          ER = 0.
          GAM_LHT = 1.
          GAM_AGE = 1.
          GAM_SMT = 1.
          GAM_TLD = 1.
          GAM_TLI = 1.
          CFNO    = 1.
          CFNOG   = 1.

          CALL GAMMA_LAI (NCOLS, NROWS, LAIc, GAM_LHT)

          CALL GAMMA_A (IDATE, ITIME, NCOLS, NROWS,  VNAME3D(S),   
     &                  LAIp, LAIc, TSTLEN, D_TEMP, GAM_AGE)

          CALL GAMMA_S (NCOLS, NROWS, GAM_SMT)

          WRITE(MESG,1030) 'Entering CANOPY: ',S,IDATE,ITIME
          CALL M3MESG( MESG )
          DO I=1, NCOLS
          DO J=1,NROWS
            ADJUST_FACTOR_LD = 0.0
            ADJUST_FACTOR_LI = 0.0
            TotalPFT = 0.0
            DO I_PFT = 1,NRTYP   !canopy types
              TotalPFT = TotalPFT + PFTF(I_PFT,I,J) * 0.01
            ENDDO   ! ENDDO I_PFT

            IF (TotalPFT .GT. 0.0) THEN

              DO I_PFT = 1,NRTYP   !canopy types
                IF (PFTF(I_PFT,I,J) .NE. 0.0) THEN
                CALL GAMME_CE(IDATE,ITIME,LAT(I,J),LONG(I,J),
     &                 TEMP(I,J),D_TEMP(I,J), D_TEMP(I,J),    
     &                 PPFD(I,J),D_PPFD(I,J), D_PPFD(I,J),   
     &                 WIND(I,J),QV(I,J),             
     &                 I_PFT,LAIc(i,j),PRES(I,J),DI(I,J),
     &                 NrCha,NrTyp,Canopychar, VNAME3D(s),
     &                 GAMMA_TD,GAMMA_TI)
                ENDIF
                ADJUST_FACTOR_LD = ADJUST_FACTOR_LD +     
     &                             0.01*PFTF(I_PFT,I,J)*GAMMA_TD
                ADJUST_FACTOR_LI = ADJUST_FACTOR_LI +    
     &                             0.01*PFTF(I_PFT,I,J)*GAMMA_TI

              ENDDO   ! ENDDO I_PFT
              GAM_TLD(I,J) = ADJUST_FACTOR_LD/TotalPFT
              GAM_TLI(I,J) = ADJUST_FACTOR_LI/TotalPFT

            ELSEIF (TotalPFT .EQ. 0.0) THEN

              GAM_TLD(I,J) = 1.0
              GAM_TLI(I,J) = 1.0

            ELSEIF (TotalPFT .LT. 0.0) THEN

              CALL M3ERR(PROGNAME,IDATE,ITIME,
     &              'TotalPFT is less than 0.0',.TRUE.)

            ENDIF

          ENDDO   ! ENDDO J
          ENDDO   ! ENDDO I
          WRITE(MESG,1030) 'Exited CANOPY: ',S,IDATE,ITIME
          CALL M3MESG( MESG )

          NMAP = INDEX1( VNAME3D(s),N_LDF_SPC,LDF_SPC)
          LDF  = LDF_FCT(NMAP)
          NMAP = INDEX1( VNAME3D(s),N_RHO_SPC,RHO_SPC)
          RHO  = RHO_FCT(NMAP)

!...  Calculate emission
          ER(:,:) = GAM_AGE * GAM_SMT * RHO *      
     &                ((1-LDF) * GAM_TLI * GAM_LHT +   
     &                LDF * GAM_TLD)
          WHERE(ER(:,:).GT.0.0) 
           NON_DIMGARMA (S,:,:) = ER(:,:)
          ELSEWHERE
           NON_DIMGARMA (S,:,:) = 0.0
          ENDWHERE
        ENDDO

!...  Estimate CFNO and CFNOG
        WRITE(MESG,1030) 'Estimating soil NOx adj: ',T,IDATE,ITIME
        CALL M3MESG( MESG )
        IF ( READ3(MGNMET,'SOIM1', ALLAYS3,IDATE,ITIME,SOILM) .AND.
     &       READ3(MGNMET,'SOIT1', ALLAYS3,IDATE,ITIME,SOILT) .AND.
     &       READ3(MGNMET,'SLTYP', ALLAYS3,IDATE,ITIME,RSTYP) ) THEN

          MESG = 'Using SOIL parameters in NOx adjustment'
          CALL M3MESG( MESG )
          LSOIL = .TRUE.
          SLTYP = INT(RSTYP)
        ELSE
          MESG = 'SOIL parameters are not available'
          CALL M3MESG( MESG )
          LSOIL = .FALSE.
        ENDIF
        CALL SOILNOX(IDATE,ITIME,NCOLS,NROWS,
     &               TEMP,LSOIL,SLTYP, SOILM, SOILT,
     &               LAIc, LAT,
     &               PRECADJ,
     &               CFNO, CFNOG )
        WRITE(MESG,1030) 'Finished soil NOx adj: ',T,IDATE,ITIME
        CALL M3MESG( MESG )

!-----------------------------------------------------------------------
!.....3) Write out the calculated ER and met data
!-----------------------------------------------------------------------
!... Write emission to file
        WRITE(MESG,1030) 'Writing emission at ',T,IDATE,ITIME
        CALL M3MESG( MESG )
        DO S = 1,NEMIS
          IF (.NOT. WRITE3(MGNERS, VNAME3D(s),IDATE,ITIME,
     &              NON_DIMGARMA (S,:,:))) THEN
            CALL NAMEVAL (MGNERS, MESG)  ! get input file name and path
            MESG = 'Error writing to file: '//TRIM(MESG)
            CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
          ENDIF

        ENDDO ! End loop for emission species (S)

! #1
        IF ( .NOT. WRITE3(MGNERS,VNAME3D(NEMIS+1),IDATE,ITIME,
     &                D_TEMP(:,:))) THEN
          CALL NAMEVAL (MGNERS, MESG)  ! get input file name and path
          MESG = 'Error writing to file: '//TRIM(MESG)
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

! #2
        IF ( .NOT. WRITE3(MGNERS,VNAME3D(NEMIS+2),IDATE,ITIME,
     &               D_PPFD(:,:))) THEN
          CALL NAMEVAL (MGNERS, MESG)  ! get input file name and path
          MESG = 'Error writing to file: '//TRIM(MESG)
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

! #3
        IF ( .NOT. WRITE3(MGNERS,VNAME3D(NEMIS+3),IDATE,ITIME,
     &              LAT(:,:))) THEN
          CALL NAMEVAL (MGNERS, MESG)  ! get input file name and path
          MESG = 'Error writing to file: '//TRIM(MESG)
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
         ENDIF

! #4
        IF ( .NOT. WRITE3(MGNERS,VNAME3D(NEMIS+4),IDATE,ITIME,
     &               LONG(:,:))) THEN
          CALL NAMEVAL (MGNERS, MESG)  ! get input file name and path
          MESG = 'Error writing to file: '//TRIM(MESG)
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

! #5
        IF ( .NOT. WRITE3(MGNERS,VNAME3D(NEMIS+5),IDATE,ITIME,
     &               CFNO(:,:))) THEN
          CALL NAMEVAL (MGNERS, MESG)  ! get input file name and path
          MESG = 'Error writing to file: '//TRIM(MESG)
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

! #6
        IF ( .NOT. WRITE3(MGNERS,VNAME3D(NEMIS+6),IDATE,ITIME,
     &               CFNOG(:,:))) THEN
          CALL NAMEVAL (MGNERS, MESG)  ! get input file name and path
          MESG = 'Error writing to file: '//TRIM(MESG)
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

! #7
        IF ( .NOT. WRITE3(MGNERS,VNAME3D(NEMIS+7),IDATE,ITIME,
     &               SLTYP(:,:))) THEN
          CALL NAMEVAL (MGNERS, MESG)  ! get input file name and path
          MESG = 'Error writing to file: '//TRIM(MESG)
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

        CALL NEXTIME( IDATE, ITIME, TSTEP )
      ENDDO ! End loop for time step (T)

!... Exit and close file
      CALL M3EXIT(PROGNAME,0,0,' ',0)

      DEALLOCATE ( ER      )   ! output emission buffer
      DEALLOCATE ( NON_DIMGARMA)
      DEALLOCATE ( LAT     )   ! input latitude of grid cell
      DEALLOCATE ( LONG    )   ! input longitude of grid cell

      DEALLOCATE ( LAIp    )   ! previous monthly LAI
      DEALLOCATE ( LAIc    )   ! current monthly LAI

      DEALLOCATE ( TEMP    )   ! input hourly temperature (K)
      DEALLOCATE ( PPFD    )   ! calculated PAR (umol/m2.s)

      DEALLOCATE ( D_PPFD  )   ! daily PAR (umol/m2.s)
      DEALLOCATE ( D_TEMP  )   ! input daily temperature (K)

      DEALLOCATE ( GAM_LHT )   ! LAI correction factor
      DEALLOCATE ( GAM_AGE )   ! leaf age correction factor
      DEALLOCATE ( GAM_SMT )   ! Soil moilture correction factor

      DEALLOCATE ( WIND    )   
      DEALLOCATE ( PRES    )   
      DEALLOCATE ( DI      )   
      DEALLOCATE ( SLTYP   )
      DEALLOCATE ( RSTYP   )
      DEALLOCATE ( SOILM   )
      DEALLOCATE ( SOILT   )
      DEALLOCATE ( QV      )
      DEALLOCATE ( PRECADJ )
      DEALLOCATE ( PFTF    )
      DEALLOCATE ( GAM_TLD )
      DEALLOCATE ( GAM_TLI )
      DEALLOCATE ( CFNO    )
      DEALLOCATE ( CFNOG   )
!--=====================================================================
!...  FORMAT
!--=====================================================================
1000  FORMAT( A )
1010  FORMAT( 43( A, :, I8, :, 1X ) )
1020  FORMAT (A40,I8,X,I8,X,I8)
1030  FORMAT (A20,I8,X,I8,X,I8)

!--=====================================================================
!...  End program
!--=====================================================================
      END PROGRAM EMPROC
