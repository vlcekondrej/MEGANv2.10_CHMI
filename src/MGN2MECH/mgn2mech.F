      PROGRAM MGN2MECH

!***********************************************************************
!   This program does chemical speciation and MECHANISM conversion.
!   The output from megan.f is converted from 20 to 150 species which
!   are then lumped according to the MECHANISM assigned in the run script.  
!   The program loops through all timesteps of the input file.
!
!   Procedure
!   1) File set up and assign I/O parameters
!   2) Conversion from MGN 20 to speciated 150
!   3) Conversion from speciated species to MECHANISM species
!   4) Convert to tonne/hour if needed
!
!   The input file gives variables in units of g-species/sec.
!   All outputs are in mole/sec or tonne/hr depending on assignment.
!
!
!   INPUT:
!           1) MEGAN output (netCDF-ioapi)
!
!   OUTPUT:
!           1) MEGAN speciation or MECHANISM species (netCDF-ioapi)
!
!   Requirement:
!      Requires libnetcdf.a and libioapi.a to compile
!
!      setenv MGERFILE    <DEFANGED_input MEGAN output for emission activity factors>  
!      setenv OUTPFILE    <output speciated emission>
!
!   CALLS:  CHECKMEM
!
!   Originally created by Jack Chen 11/04 for MEGAN v.0
!   For MEGAN v2.0 created by Tan 12/01/06
!   For MEGAN v2.1 created by Xuemei Wang 11/04/07
!   For MEGAN v2.1 to use 150 species created by Xuemei Wang 09/30/09
!
!  History:
!  08/14/07 Tan    - Move to MEGANv2.02 with no update
!  08/29/07 modified by A. Guenther to correct error in assigning   
!           emission factor. This version is called MEGANv2.03
!  10/29/07 modified by A. Guenther to correct omission of diurnal variation    
!           factor. This version is called MEGANv2.04
!  11/04/07 modified by Xuemei Wang to give two options for MAP or lookup table for
!           the emission factors. Also gives options for different chemical MECHANISMs
!           in the code: user modifies the external script to assign MECHANISM.
!           This version is called MEGANv2.1.0
!  06/04/08 modified by J. Lee-Taylor to accept vegetation-dependent speciation factors
!           in table format (RESHAPE tables) rather than from DATA statements.
!  09/30/08  modified by Xuemei Wang to give options for input file and test different mechanisms
!  09/27/11  Tan&Xuemei MEGANv2.10 includes soil NOx adjustment and a lot of updates
!***********************************************************************

      IMPLICIT NONE

!...  INCLUDES:
      INCLUDE 'PARMS3.EXT'        !  I/O API parameters
      INCLUDE 'IODECL3.EXT'       !  I/O API function declarations
      INCLUDE 'FDESC3.EXT'        !  I/O API file description data structures.

      INCLUDE 'SPC_NOCONVER.EXT'
      INCLUDE 'EFS_PFT.EXT'
      INCLUDE 'MAP_CV2MECHANISM.EXT'     !  Conversion map and factors for different MECHANISM 
      INCLUDE 'PFT_MGN.EXT'
      INCLUDE 'SPC_SAPRC99.EXT'
      INCLUDE 'MAP_CV2SAPRC99.EXT'
      INCLUDE 'SPC_SAPRC99X.EXT'
      INCLUDE 'MAP_CV2SAPRC99X.EXT'
      INCLUDE 'SPC_SAPRC99Q.EXT'
      INCLUDE 'MAP_CV2SAPRC99Q.EXT'
      INCLUDE 'SPC_CB05.EXT'
      INCLUDE 'MAP_CV2CB05.EXT'
      INCLUDE 'SPC_CB6.EXT'
      INCLUDE 'MAP_CV2CB6.EXT'
      INCLUDE 'SPC_SOAX.EXT'
      INCLUDE 'MAP_CV2SOAX.EXT'

C...........   EXTERNAL FUNCTIONS and their descriptions:
      INTEGER      PROMPTFFILE, STR2INT
      REAL         STR2REAL
      INTEGER      INDEX1
      EXTERNAL     PROMPTFFILE,STR2INT
      EXTERNAL     STR2REAL
      EXTERNAL     INDEX1
      LOGICAL      ENVYN
      EXTERNAL     ENVYN
      LOGICAL      DSCGRID
      EXTERNAL     DSCGRID


!...  Program I/O files
! Program name
      CHARACTER*16 :: PROGNAME = 'MGN2MECH'
! Input MEGAN ER file
      CHARACTER*16 :: MGNERS   = 'MGNERS'    ! Input MEGAN ER file logical name
! Netcdf file
      CHARACTER*16 :: EFMAPS   = 'EFMAPS'    !  EFMAP input file  name
      CHARACTER*16 :: PFTS16   = 'PFTS16'    ! Input PFT file logical
! Output file
      CHARACTER*16 :: MGNOUT   = 'MGNOUT'    ! Output file logical name
! Parameters for file units
      INTEGER :: logdev                      ! Logfile unit number

!...  Program I/O parameters
!...  External parameters

! from run script
      LOGICAL :: TONPHR      = .FALSE.  ! Output in tons/hr flag
      LOGICAL :: CONVERSION  = .FALSE.  ! Mechanism conversion flag
      !LOGICAL :: NONPFT_EF   = .FALSE.

      CHARACTER*16 :: MECHANISM              ! Mechanism name

! I/O API file parameters
      INTEGER       jdate          ! Looping date YYYYDDD
      INTEGER       jtime          ! Looping time HHMMSS
      INTEGER       ncols          ! Number of columns
      INTEGER       nrows          ! Number of rows
      INTEGER       mxrec          ! Total number of time steps
      INTEGER       sdate          ! Start date YYYYDDD
      INTEGER       stime          ! Start time HHMMSS
      INTEGER       tstep          ! Time step


!...  Internal parameters
! internal paramters (status and buffer)
      INTEGER   ios                          ! i/o status
      CHARACTER*256 :: MESG                  ! Message buffer

      INTEGER :: cid, inx, iny               ! Input grid x and y
      INTEGER :: MXPFT

! local variables and their descriptions:
      CHARACTER*16  :: GDNAM
      CHARACTER*16  :: CNAME        ! Coord name

      INTEGER :: t, s, I, N, C, R             ! Counters
      INTEGER :: nmpmg, nmpsp, nmpmc          ! Counters
      REAL, ALLOCATABLE :: LAT( :,: )         ! Latitude of grid cell
      REAL, ALLOCATABLE :: inper(:,:,:)       ! Input emission buffer
      REAL, ALLOCATABLE :: tmper(:,:,:)       ! Temp emission buffer
      REAL, ALLOCATABLE :: outer(:,:,:)       ! Output emission buffer
      REAL, ALLOCATABLE :: EF(:,:,:)          ! Emission factor  
      REAL, ALLOCATABLE :: TMP1(:,:),TMP2(:,:),TMP3(:,:)
      REAL              :: TMO1,TMO2,TMO3
      REAL, ALLOCATABLE :: CFNOG(:,:), CFNO(:,:)
      REAL, ALLOCATABLE :: PFT(:,:,:)         ! PFT factor array
      INTEGER :: n_scon_spc 
      INTEGER :: NVAR
      INTEGER       GDAY, GLEN, INO, IISOP
      REAL          garea , G2CRatio
      INTEGER, ALLOCATABLE:: GARR(:,:)

      INTEGER, ALLOCATABLE ::  spmh_map(:),mech_map(:)   ! speciated species name
      CHARACTER*16, ALLOCATABLE :: mech_spc(:)
      REAL, ALLOCATABLE :: conv_fac(:)        
      REAL,ALLOCATABLE :: mech_mwt(:)
      
      REAL, PARAMETER :: ug2g = 1E-6          ! convert microgram to metric gram
      REAL, PARAMETER :: g2tonne = 1E-6       ! convert microgram to metric ton
      REAL, PARAMETER :: hr2sec = 3600        ! convert hr to second
      REAL, PARAMETER :: n2no   = 2.142857    ! convert hr to second
      

C***********************************************************************

!=======================================================================
!...  Begin program
!=======================================================================


!-----------------------------------------------------------------------
!.....1) File set up and assign I/O parameters
!-----------------------------------------------------------------------
!...  Initialize log file unit
      logdev = init3()

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

      IF ( .NOT. OPEN3( MGNERS, FSREAD3, PROGNAME ) ) THEN
          MESG = 'Could not open input file '//TRIM(MGNERS)
          CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      ! Check grid
      IF ( .NOT. FILCHK3 ( MGNERS,
     &              GRDDED3, NCOLS3D, NROWS3D, 1, NTHIK3D))  THEN
         MESG = 'MGNERS has differenet grid definition'
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      IF ( .NOT. DESC3( MGNERS ) ) THEN
          MESG = 'Could not get description of '//TRIM(MGNERS)
          CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF

      sdate = sdate3d
      stime = stime3d
      tstep = tstep3d
      mxrec = mxrec3d
      ncols = ncols3d
      nrows = nrows3d
      garea = xcell3d * ycell3d
      jdate = sdate
      jtime = stime
      print *,' Area per grid cell (m2): ',garea

!...  Set output MECHANISM variable variables
! Check to run MECHANISM conversion
      MESG = 'Convert to any ATM MECHANISM'
      CONVERSION = ENVYN ( 'RUN_CONVERSION', MESG, .FALSE., ios )

! Check tonne per hour output
      MESG = 'Output in Tonne per hour'
      TONPHR = ENVYN ( 'SPCTONHR', MESG, .FALSE., ios )
      !MESG = 'Emission factor with based on PFT or not'
      !NONPFT_EF = ENVYN ( 'SPC_EF', MESG, .FALSE., IOS )


! Set attribute and variables for output
      IF ( CONVERSION ) THEN
        MESG = 'Mechanism Name'
        CALL ENVSTR( 'MECHANISM', MESG, 'SAPRCII', MECHANISM, ios )

        SELECT CASE ( TRIM(MECHANISM) )
          CASE ('SAPRCII')
            n_scon_spc = n_saprcii
            NVARS3D = n_saprcii_spc
          CASE ('SAPRC99')
            n_scon_spc = n_saprc99
            NVARS3D = n_saprc99_spc
          CASE ('RADM2')
            n_scon_spc = n_radm2
            NVARS3D = n_radm2_spc
          CASE ('RACM')
            n_scon_spc = n_racm
            NVARS3D = n_racm_spc
          CASE ('CBMZ')
            n_scon_spc = n_cbmz
            NVARS3D = n_cbmz_spc
          CASE ('SAPRC99X')
            n_scon_spc = n_saprc99_x
            NVARS3D = n_saprc99_x_spc
          CASE ('SAPRC99Q')
            n_scon_spc = n_saprc99_q
            NVARS3D = n_saprc99_q_spc
          CASE ('CB05')
            n_scon_spc = n_cb05
            NVARS3D = n_cb05_spc
          CASE ('CB6')
            n_scon_spc = n_cb6
            NVARS3D = n_cb6_spc
          CASE ('SOAX')
            n_scon_spc = n_soax
            NVARS3D = n_soax_spc
          CASE DEFAULT
            MESG = 'Error: Mechanism conversion, invalid MECHANISM: '
     &            //TRIM(MECHANISM)
            CALL M3EXIT(PROGNAME, 0, 0,MESG,2)
        ENDSELECT
 
        ALLOCATE ( spmh_map(n_scon_spc), STAT = ios )
          CALL CHECKMEM ( ios, 'spmh_map', PROGNAME )
        ALLOCATE ( mech_map(n_scon_spc), STAT = ios )
          CALL CHECKMEM ( ios, 'mech_map', PROGNAME )
        ALLOCATE ( conv_fac(n_scon_spc), STAT = ios )
          CALL CHECKMEM ( ios, 'conv_fac', PROGNAME )
        ALLOCATE ( mech_spc(NVARS3D ), STAT = ios )
          CALL CHECKMEM ( ios, 'mech_spc', PROGNAME )
        ALLOCATE ( mech_mwt(NVARS3D ), STAT = ios )
          CALL CHECKMEM ( ios, 'mech_mwt', PROGNAME )

        SELECT CASE ( TRIM(MECHANISM) )

          CASE ('SAPRCII')
            spmh_map(1:n_scon_spc) = spmh_map_saprcii(1:n_scon_spc)
            mech_map(1:n_scon_spc) = mech_map_saprcii(1:n_scon_spc)
            conv_fac(1:n_scon_spc) = conv_fac_saprcii(1:n_scon_spc)
            mech_spc(1:NVARS3D)    = mech_spc_saprcii(1:NVARS3D)
            mech_mwt(1:NVARS3D)    = mech_mwt_saprcii(1:NVARS3D)
            VNAME3D(1:NVARS3D)     = mech_spc(1:NVARS3D)
          CASE ('SAPRC99')
            spmh_map(1:n_scon_spc) = spmh_map_saprc99(1:n_scon_spc)
            mech_map(1:n_scon_spc) = mech_map_saprc99(1:n_scon_spc)
            conv_fac(1:n_scon_spc) = conv_fac_saprc99(1:n_scon_spc)
            mech_spc(1:NVARS3D)    = mech_spc_saprc99(1:NVARS3D)
            mech_mwt(1:NVARS3D)    = mech_mwt_saprc99(1:NVARS3D)
            VNAME3D(1:NVARS3D)     = mech_spc(1:NVARS3D)
          CASE ('RADM2')
            spmh_map(1:n_scon_spc) = spmh_map_radm2(1:n_scon_spc)
            mech_map(1:n_scon_spc) = mech_map_radm2(1:n_scon_spc)
            conv_fac(1:n_scon_spc) = conv_fac_radm2(1:n_scon_spc)
            mech_spc(1:NVARS3D)    = mech_spc_radm2(1:NVARS3D)
            mech_mwt(1:NVARS3D)    = mech_mwt_radm2(1:NVARS3D)
            VNAME3D(1:NVARS3D)     = mech_spc(1:NVARS3D)
          CASE ('RACM')
            spmh_map(1:n_scon_spc) = spmh_map_racm(1:n_scon_spc)
            mech_map(1:n_scon_spc) = mech_map_racm(1:n_scon_spc)
            conv_fac(1:n_scon_spc) = conv_fac_racm(1:n_scon_spc)
            mech_spc(1:NVARS3D)    = mech_spc_racm(1:NVARS3D)
            mech_mwt(1:NVARS3D)    = mech_mwt_racm(1:NVARS3D)
            VNAME3D(1:NVARS3D)     = mech_spc(1:NVARS3D)
          CASE ('CBMZ')
            spmh_map(1:n_scon_spc) = spmh_map_cbmz(1:n_scon_spc)
            mech_map(1:n_scon_spc) = mech_map_cbmz(1:n_scon_spc)
            conv_fac(1:n_scon_spc) = conv_fac_cbmz(1:n_scon_spc)
            mech_spc(1:NVARS3D)    = mech_spc_cbmz(1:NVARS3D)
            mech_mwt(1:NVARS3D)    = mech_mwt_cbmz(1:NVARS3D)
            VNAME3D(1:NVARS3D)     = mech_spc(1:NVARS3D)
          CASE ('SAPRC99X')
            spmh_map(1:n_scon_spc) = spmh_map_saprc99_X(1:n_scon_spc)
            mech_map(1:n_scon_spc) = mech_map_saprc99_X(1:n_scon_spc)
            conv_fac(1:n_scon_spc) = conv_fac_saprc99_X(1:n_scon_spc)
            mech_spc(1:NVARS3D)    = mech_spc_saprc99_X(1:NVARS3D)
            mech_mwt(1:NVARS3D)    = mech_mwt_saprc99_X(1:NVARS3D)
            VNAME3D(1:NVARS3D)     = mech_spc(1:NVARS3D)
          CASE ('SAPRC99Q')
            spmh_map(1:n_scon_spc) = spmh_map_saprc99_Q(1:n_scon_spc)
            mech_map(1:n_scon_spc) = mech_map_saprc99_Q(1:n_scon_spc)
            conv_fac(1:n_scon_spc) = conv_fac_saprc99_Q(1:n_scon_spc)
            mech_spc(1:NVARS3D)    = mech_spc_saprc99_Q(1:NVARS3D)
            mech_mwt(1:NVARS3D)    = mech_mwt_saprc99_Q(1:NVARS3D)
            VNAME3D(1:NVARS3D)     = mech_spc(1:NVARS3D)
          CASE ('CB05')
            spmh_map(1:n_scon_spc) = spmh_map_cb05(1:n_scon_spc)
            mech_map(1:n_scon_spc) = mech_map_cb05(1:n_scon_spc)
            conv_fac(1:n_scon_spc) = conv_fac_cb05(1:n_scon_spc)
            mech_spc(1:NVARS3D)    = mech_spc_cb05(1:NVARS3D)
            mech_mwt(1:NVARS3D)    = mech_mwt_cb05(1:NVARS3D)
            VNAME3D(1:NVARS3D)     = mech_spc(1:NVARS3D)
          CASE ('CB6')
            spmh_map(1:n_scon_spc) = spmh_map_cb6(1:n_scon_spc)
            mech_map(1:n_scon_spc) = mech_map_cb6(1:n_scon_spc)
            conv_fac(1:n_scon_spc) = conv_fac_cb6(1:n_scon_spc)
            mech_spc(1:NVARS3D)    = mech_spc_cb6(1:NVARS3D)
            VNAME3D(1:NVARS3D)     = mech_spc(1:NVARS3D)
          CASE ('SOAX')
            spmh_map(1:n_scon_spc) = spmh_map_soax(1:n_scon_spc)
            mech_map(1:n_scon_spc) = mech_map_soax(1:n_scon_spc)
            conv_fac(1:n_scon_spc) = conv_fac_soax(1:n_scon_spc)
            mech_spc(1:NVARS3D)    = mech_spc_soax(1:NVARS3D)
            mech_mwt(1:NVARS3D)    = mech_mwt_soax(1:NVARS3D)
            VNAME3D(1:NVARS3D)     = mech_spc(1:NVARS3D)
        ENDSELECT
      ELSEIF ( .NOT. CONVERSION ) THEN
        print*,'No conversion'
        NVARS3D = n_spca_spc
        VNAME3D(1:NVARS3D) = spca_spc(1:NVARS3D)
      ELSE
         MESG = 'Error: Conversion flag is not assigned.'
         CALL M3EXIT(PROGNAME, 0, 0,MESG,2)
      ENDIF

! Change the unit according to TONPHR flag
      IF ( TONPHR ) THEN
         UNITS3D(1:NVARS3D) = 'tons/hr'
      ELSE
         UNITS3D(1:NVARS3D) = 'mol/s'
      ENDIF
      VDESC3D(1:NVARS3D) = ' '
      VTYPE3D(1:NVARS3D) = m3real
      DO s = 1, NVARS3D
         PRINT*,'Output variable:',VNAME3D(s),UNITS3D(s)
      ENDDO
      ftype3d = grdded3
      sdate3d = sdate
      stime3d = stime
      tstep3d = tstep
      mxrec3d = mxrec

      IF ( CONVERSION ) THEN
      MESG = 'Mechanism Name'
      CALL ENVSTR( 'MECHANISM', MESG, 'SAPRCII', MECHANISM, ios )
         SELECT CASE ( TRIM(MECHANISM) )
         CASE ('SAPRCII')
            FDESC3D( 1 ) = 'Chemical specation of MEGAN output to '
     &                     //'SAPRCII Species:'
         CASE ('SAPRC99')
            FDESC3D( 1 ) = 'Chemical specation of MEGAN output to '
     &                     //'SAPRC99 Species:'
         CASE ('RADM2')
            FDESC3D( 1 ) = 'Chemical specation of MEGAN output to '
     &                     //'RADM2 Species:'
         CASE ('RACM')
            FDESC3D( 1 ) = 'Chemical specation of MEGAN output to '
     &                     //'RACM Species:'
         CASE ('CBMZ')
            FDESC3D( 1 ) = 'Chemical specation of MEGAN output to '
     &                     //'CBMZ Species:'
         CASE ('SAPRC99X')
            FDESC3D( 1 ) = 'Chemical specation of MEGAN output to '
     &                     //'SAPRC99X Species:'
         CASE ('SAPRC99Q')
            FDESC3D( 1 ) = 'Chemical specation of MEGAN output to '
     &                     //'SAPRC99Q Species:'
         CASE ('CB05')
            FDESC3D( 1 ) = 'Chemical specation of MEGAN output to '
     &                     //'CB05 Species:'
         CASE ('CB6')
            FDESC3D( 1 ) = 'Chemical specation of MEGAN output to '
     &                     //'CB6 Species:'
         CASE ('SOAX')
            FDESC3D( 1 ) = 'Chemical specation of MEGAN output to '
     &                     //'SOAX Species:'
         CASE DEFAULT
            MESG = 'Error: Mechanism conversion, invalid MECHANISM: '
     &             //TRIM(MECHANISM)
            CALL M3EXIT(PROGNAME,JDATE,JTIME,MESG,2)
         ENDSELECT

      ELSEIF ( .NOT. CONVERSION ) THEN
         FDESC3D( 1 ) = 'Chemical specation of MEGAN output to '
     &                  //'MEGAN 150 species'
      ELSE
         MESG = 'Error: Conversion flag is not assigned.'
         CALL M3EXIT(PROGNAME, 0, 0,MESG,2)

      ENDIF

      CALL NAMEVAL ( MGNERS , MESG )  ! get input file name and path
      FDESC3D( 2 ) = 'Input MEGAN file: '//TRIM(MESG)

      VNAME3D(NVARS3D+1) = 'GDAY'
      UNITS3D(NVARS3D+1) = 'day'
      VDESC3D(NVARS3D+1) = ' '
      VTYPE3D(NVARS3D+1) = m3int
      NVAR = NVARS3D
      NVARS3D = NVARS3D + 1
      IF ( .NOT. OPEN3( MGNOUT, FSCREA3, PROGNAME ) ) THEN
        MESG = 'Could not open file '//TRIM(MGNOUT)
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF

!...  Allocate memory
      ALLOCATE ( EF   ( N_MGN_SPC, ncols, nrows ), STAT = ios )
        CALL CHECKMEM ( ios, 'EF',    PROGNAME )

      ALLOCATE ( TMP1 (  ncols, nrows ), STAT = ios )
        CALL CHECKMEM ( ios, 'tmp1',  PROGNAME )

      ALLOCATE ( TMP2 (  ncols, nrows ), STAT = ios )
        CALL CHECKMEM ( ios, 'tmp2',  PROGNAME )

      ALLOCATE ( TMP3 (  ncols, nrows ), STAT = ios )
        CALL CHECKMEM ( ios, 'tmp3',  PROGNAME )

      ALLOCATE ( LAT  (  ncols, nrows ), STAT = ios )
        CALL CHECKMEM ( ios, 'LAT',   PROGNAME )

      ALLOCATE ( CFNO (  ncols, nrows ), STAT = ios )
        CALL CHECKMEM ( ios, 'CFNO',  PROGNAME )

      ALLOCATE ( CFNOG(  ncols, nrows ), STAT = ios )
        CALL CHECKMEM ( ios, 'CFNOG', PROGNAME )

      ALLOCATE ( PFT  ( MXPFT, ncols, nrows ), STAT = ios )
        CALL CHECKMEM ( ios, 'PFT',   PROGNAME )

      ALLOCATE ( GARR ( ncols, nrows ), STAT = ios )
        CALL CHECKMEM ( ios, 'GARR', PROGNAME )


!...  Read EFMAPS

      IF ( .NOT. READ3(EFMAPS,'LAT',1,1,1,LAT)) THEN
         MESG = 'Error reading LAT'
         CALL M3EXIT(PROGNAME,jdate,jtime,MESG,2)
      ENDIF

      DO N = 1, N_MGN_SPC
        IF ( .NOT. READ3(EFMAPS,'EF_'//MGN_SPC(N),1,1,1,EF(N,:,:))) THEN
           MESG = 'Error reading '//TRIM('EF_'//MGN_SPC(N))
           CALL M3EXIT(PROGNAME,jdate,jtime,MESG,2)
        ENDIF
      ENDDO

!...  Read PFTS16
      DO N = 1, MXPFT
         IF ( .NOT. READ3(PFTS16,'PFTS',1,0,(N-1)*10000,
     &      PFT(N,:,:))) THEN
            MESG = 'Error reading PFTS'
            CALL M3EXIT(PROGNAME,0,(N-1)*10000,MESG,2)
         ENDIF
      ENDDO
       
!.....2) Conversion from MGN 20 to speciated 150
!-----------------------------------------------------------------------
!...  Allocate memory
      ALLOCATE ( inper(N_MGN_SPC, ncols, nrows ), STAT = ios )
        CALL CHECKMEM ( ios, 'inper', PROGNAME )
      ALLOCATE ( tmper( n_spca_spc, ncols, nrows ), STAT = ios )
        CALL CHECKMEM ( ios, 'tmper', PROGNAME )
      ALLOCATE ( outer( NVAR, ncols, nrows ), STAT = ios )
        CALL CHECKMEM ( ios, 'outer', PROGNAME )


      INO = INDEX1('NO',N_MGN_SPC,MGN_SPC)
      IISOP = INDEX1('ISOP',N_MGN_SPC,MGN_SPC)
!...  Loop through time
      jdate = sdate
      jtime = stime
      DO t = 1, mxrec
        inper = 0.
        tmper = 0.
        outer = 0.

        DO N = 1, N_MGN_SPC
          IF ( .NOT. READ3(MGNERS,TRIM(MGN_SPC(N)),1,
     &                      jdate,jtime,inper(N,:,:))) THEN
            MESG = 'Error reading '//TRIM(MGN_SPC(nmpmg))
            CALL M3EXIT(PROGNAME,jdate,jtime,MESG,2)
          ENDIF
        ENDDO

        IF ( .NOT. READ3(MGNERS,'CFNO',1,jdate,jtime,CFNO(:,:))) THEN
          MESG = 'Error reading CFNO'
          CALL M3EXIT(PROGNAME,jdate,jtime,MESG,2)
        ENDIF

        IF( .NOT. READ3(MGNERS,'CFNOG',1,jdate,jtime,CFNOG(:,:))) THEN
          MESG = 'Error reading CFNOG'
          CALL M3EXIT(PROGNAME,jdate,jtime,MESG,2)
        ENDIF

        DO s = 1, N_SMAP_SPC
            
          nmpmg = mg20_map(s)
          nmpsp = spca_map(s)
          IF (t .EQ. 1) THEN
            MESG='Convert '//MGN_SPC(NMPMG)//' to '//SPCA_SPC(NMPSP)
            CALL M3MESG( MESG )
          ENDIF

         IF ( nmpmg .NE. INO ) then
           !...  Not NO
           IF ( ef_all(1,nmpmg) .LT. 0.0) THEN
             !... Use EFMAPS
             MESG = 'Use EFMAPS for '//MGN_SPC(NMPMG)
             CALL M3MESG( MESG )
             TMP1(:,:) = 0.0
             TMP2(:,:) = 0.0
             DO I = 1,N_MGN_PFT
               TMP1 = TMP1 + PFT(I,:,:)
               TMP2 = TMP2 + effs_all(I,nmpsp) * PFT(I,:,:)
             ENDDO
             WHERE(TMP1(:,:) .EQ. 0.0)
               tmper(nmpsp,:,:) = 0.0
             ELSEWHERE
               tmper(nmpsp,:,:) = inper(nmpmg,:,:) * EF(nmpmg,:,:)
     &                             * TMP2(:,:) / TMP1(:,:)
             ENDWHERE
           ELSE
             !... Use PFT-EF
             TMP3(:,:) = 0.0
             DO I = 1,N_MGN_PFT
               TMP3 = TMP3 + ef_all(I,nmpmg) * effs_all(I,nmpsp)
     &                         * PFT(I,:,:)/100.0
             ENDDO
             tmper(nmpsp,:,:) = inper(nmpmg,:,:) * TMP3(:,:)
           ENDIF
         ELSE IF ( nmpmg .EQ. INO ) then
!!-----------------NO Stuff-----------------------
           !...  NO when nmpmg = 17
           MESG = 'Special loop for '//MGN_SPC(NMPMG)
           CALL M3MESG( MESG )
           IF (ef_all(1,INO) .LT. 0.0) THEN
             MESG = 'Use EFMAPS for '//MGN_SPC(NMPMG)
             CALL M3MESG( MESG )
           ENDIF
           DO C = 1,ncols
           DO R = 1,nrows
           CALL GROWSEASON(jdate,LAT(C,R),GDAY,GLEN)
           GARR(C,R) = GDAY
           ! check for growing season
           IF (GDAY .EQ. 0) THEN
             ! non growing season
             ! CFNOG for everywhere
             ! Override crop with grass warm = 14
             IF (ef_all(1,INO) .LT. 0.0) THEN
               ! with EFMAPS
               TMO1 = 0.0
               TMO2 = 0.0
               DO I = 1,14
                 TMO1 = TMO1 + PFT(I,C,R)
                 TMO2 = TMO2 + effs_all(I,nmpsp) * PFT(I,C,R)
               ENDDO
               DO I = 15,N_MGN_PFT
                 TMO1 = TMO1 + PFT(I,C,R)
                 G2CRatio = ef_all(14,INO)/ef_all(I,INO)
                 TMO2 = TMO2 + effs_all(I,nmpsp) * PFT(I,C,R) * G2CRatio
               ENDDO
               IF (TMO1 .EQ. 0.0) THEN
                 tmper(nmpsp,C,R) = 0.0
               ELSE
!                 tmper(nmpsp,C,R) = inper(INO,C,R) * EF(INO,C,R) 
!     &                            * CFNOG(C,R) * TMO2 / TMO1
                 tmper(nmpsp,C,R) = inper(INO,C,R) * EF(INO,C,R)
     &                            * CFNOG(C,R) * TMO2 / TMO1
     &                            * n2no
               ENDIF
             ELSE
               ! without EFMAPS
               TMO3 = 0.0
               DO I = 1,14
                 TMO3 = TMO3 + ef_all(I,INO)*effs_all(I,nmpsp)
     &                             * PFT(I,C,R)/100.0
               ENDDO
               DO I = 15,N_MGN_PFT
                 TMO3 = TMO3 + ef_all(14,INO) * effs_all(I,nmpsp)
     &                             * PFT(I,C,R)/100.0
               ENDDO
!               tmper(nmpsp,C,R) = inper(INO,C,R) * CFNOG(C,R) * TMO3
               tmper(nmpsp,C,R) = inper(INO,C,R) * CFNOG(C,R) * TMO3
     &                            * n2no
             ENDIF
    
           ELSE IF (GDAY .GT. 0 .AND. GDAY .LE. 366) THEN
             ! growing season
             ! CFNOG for everywhere except crops
             ! CFNO for crop and corn
             IF (ef_all(1,INO) .LT. 0.0) THEN
               ! with EFMAPS
               TMO1 = 0.0
               TMO2 = 0.0
               DO I = 1,14
                 TMO1 = TMO1 + PFT(I,C,R)
                 TMO2 = TMO2 + effs_all(I,nmpsp)*PFT(I,C,R)
     &                         * CFNOG(C,R)
               ENDDO
               DO I = 15,N_MGN_PFT
                 TMO1 = TMO1 + PFT(I,C,R)
                 TMO2 = TMO2 + effs_all(I,nmpsp)*PFT(I,C,R)
     &                         *CFNO(C,R)
               ENDDO
               IF (TMO1 .EQ. 0.0) THEN
                 tmper(nmpsp,C,R) = 0.0
               ELSE
!                 tmper(nmpsp,C,R) = inper(INO,C,R) * EF(INO,C,R) *
!     &                              TMO2 / TMO1
                 tmper(nmpsp,C,R) = inper(INO,C,R) * EF(INO,C,R) *
     &                              TMO2 / TMO1
     &                              * n2no
               ENDIF
             ELSE
               ! without EFMAPS
               TMO3 = 0.0
               DO I = 1,14
                 TMO3 = TMO3 + ef_all(I,INO)*effs_all(I,nmpsp)*
     &                       PFT(I,C,R)/100.0* CFNOG(C,R)
               ENDDO
               DO I = 15,N_MGN_PFT
                 TMO3 = TMO3 + ef_all(I,INO)*effs_all(I,nmpsp)*
     &                       PFT(I,C,R)/100.0* CFNO(C,R)
               ENDDO
!               tmper(nmpsp,C,R) = inper(INO,C,R) * TMO3
               tmper(nmpsp,C,R) = inper(INO,C,R) * TMO3
     &                            * n2no
             ENDIF
           ELSE
             ! bad GDAY
             WRITE(MESG,*) 'Bad GDAY ',GDAY
             CALL M3EXIT(PROGNAME,jdate,jtime,MESG,2)
           ENDIF
           ENDDO  !DO R = 1,nrows
           ENDDO  !DO C = 1,ncols
!-----------------end of NO----------------------
         ENDIF     !IF ( nmpmg .NE. INO ) then

        ENDDO ! End species loop

!-----------------------------------------------------------------------
!.....3) Conversion from speciated species to MECHANISM species
!-----------------------------------------------------------------------
        ! convert from ug/m^2/hr to mol/m^2/hr using their MW
        DO s = 1, n_spca_spc
           tmper(s,:,:) = tmper(s,:,:) / spca_mwt(s)*ug2g
           IF (t .EQ. 1) THEN
           PRINT*, TRIM(spca_spc(s))//' = '//
     &             TRIM(spca_spc(s))//' / MW ',spca_mwt(s)
           ENDIF
        ENDDO

        IF ( CONVERSION ) THEN
          ! lumping to MECHANISM species

          DO s = 1, n_scon_spc
            nmpsp = spmh_map(s)         ! Mapping value for SPCA
            nmpmc = mech_map(s)         ! Mapping value for MECHANISM

            IF (t .EQ. 1) THEN
              MESG='Convert '//SPCA_SPC(NMPSP)//' to '//MECH_SPC(NMPMC)
              CALL M3MESG( MESG )
            ENDIF

            outer(nmpmc,:,:) = outer(nmpmc,:,:) +
     &                         (tmper(nmpsp,:,:)*conv_fac(s))*
     &                         garea/hr2sec
            ! units of these species are in mole/s
          ENDDO ! End species loop

        ELSE
          ! get all 150 species into the output array
          outer(:,:,:) = tmper(:,:,:)*
     &                         garea/hr2sec
          ! units of these species are in mole/s/grid

        ENDIF

!-----------------------------------------------------------------------
!.....4) Convert to tonne/hour if needed
!-----------------------------------------------------------------------
        IF ( TONPHR ) THEN
          IF ( CONVERSION ) THEN
            ! convert from mol/s to ton/hr
            DO s = 1, NVAR
              outer(s,:,:) = outer(s,:,:)*mech_mwt(s)*
     &                       g2tonne*hr2sec
            ENDDO
          ELSEIF ( .NOT. CONVERSION ) THEN
            ! convert from mol/s to ton/hr
            DO s = 1, NVAR
              outer(s,:,:) = outer(s,:,:)*spca_mwt(s)*
     &                       g2tonne*hr2sec
            ENDDO
          ENDIF
        ENDIF

!-----------------------------------------------------------------------
!.....5) Write out the calculated emissions
!-----------------------------------------------------------------------
        print*,'VARS=',NVAR,n_scon_spc
        DO s = 1, NVAR
          IF ( .NOT. WRITE3(MGNOUT,VNAME3D(s),jdate,jtime,
     &                                        outer(s,:,:)    )  ) THEN
             MESG = 'Error writing to file: '//TRIM(MGNOUT)
             CALL M3EXIT(PROGNAME,jdate,jtime,MESG,2)
          ENDIF
        ENDDO 

        s = NVAR + 1
        IF ( .NOT. WRITE3(MGNOUT,VNAME3D(s),jdate,jtime,
     &                                          GARR(:,:)    )  ) THEN
           MESG = 'Error writing to file: '//TRIM(MGNOUT)
           CALL M3EXIT(PROGNAME,jdate,jtime,MESG,2)
        ENDIF
!-----------------------------------------------------------------------
        CALL NEXTIME(jdate,jtime,tstep)
      ENDDO ! End time loop

!... Exit and close file
      CALL M3EXIT(PROGNAME,0,0,' ',0)

      DEALLOCATE ( spmh_map )
      DEALLOCATE ( mech_map )
      DEALLOCATE ( conv_fac )
      DEALLOCATE ( mech_spc )
      DEALLOCATE ( mech_mwt )
      DEALLOCATE (TMP1)
      DEALLOCATE (TMP2)
      DEALLOCATE (TMP3)
      DEALLOCATE (LAT)
      DEALLOCATE (EF)
      DEALLOCATE (CFNO)
      DEALLOCATE (CFNOG)
      DEALLOCATE (PFT)
      DEALLOCATE (outer)
      DEALLOCATE (inper)
      DEALLOCATE (tmper)
      DEALLOCATE (GARR)

!=======================================================================
!...  FORMAT
!=======================================================================        


!=======================================================================
!...  End program
!=======================================================================
      ENDPROGRAM
