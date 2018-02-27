
      PROGRAM TXT2IOAPI
!***********************************************************************
!  This program converts comma deliminated line input to gridded
!    netCDF-IOAPI format.
!  The grid parameter is input from GRIDDESC
!  Program always takes second and third input column as grid-x 
!    and grid-y index.
!  If you export the file from Excel to .csv file, make sure you convert 
!    the DOS text file to Unix text file with "dos2unix" command.
!
!  Requirement:
!     Requires libnetcdf.a and libioapi.a to compile
!     Also, parsline.o and chkmetem.o from SMOKE
!
!     setenv EFSTXTF <ef_map text file>
!     setenv PFTTXTF <pft    text file>
!     setenv LAITXTF <lai    text file>
!     setenv EFMAPS  <output EFS netCDF-ioapi file>
!     setenv PFTS16  <output PFT netCDF-ioapi file>
!     setenv LAIS46  <output LAI netCDF-ioapi file>
!     setenv RUN_EFS <T|F>
!     setenv RUN_LAI <T|F>
!     setenv RUN_PFT <T|F>
!     setenv GRIDDESC <grid desc file>
!
!  Call:
!     PARSLINE
!     CHECKMEM
!
!  Created by Jack Chen 11/04
!  Modified by Tan 11/13/06 for MEGAN v2.0
!
!  History:
!  08/14/07 Tan    - Move to MEGANv2.02 with no update
!  07/06/11 Tan    - Update for MEGANv2.10
!***********************************************************************

      IMPLICIT NONE

C...........   INCLUDES:
      INCLUDE 'PARMS3.EXT'    !  I/O API parameters
      INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
      INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures

      INCLUDE 'SPC_MGN.EXT'    !  MEGAN species
      INCLUDE 'LAI_MGN.EXT'    !  LAI
      INCLUDE 'PFT_MGN.EXT'    !  PFT

C...........   EXTERNAL FUNCTIONS and their descriptions:
      INTEGER      PROMPTFFILE, STR2INT
      REAL         STR2REAL
      LOGICAL      DSCGRID
      LOGICAL      ENVYN
      INTEGER      INDEX1
      EXTERNAL     PROMPTFFILE, DSCGRID
     &             STR2INT, STR2REAL, INDEX1

C.........  LOCAL VARIABLES and their descriptions:
      INTEGER       IOS            ! i/o status
      INTEGER    :: LDEV           ! unit no. for log file
      CHARACTER*256    MESG        ! message buffer

      LOGICAL       RUN_EFS
      LOGICAL       RUN_PFT
      LOGICAL       RUN_LAI

      CHARACTER*16  :: EFSTXTF = 'EFSTXTF'      ! Input EF file logical name
      CHARACTER*16  :: PFTTXTF = 'PFTTXTF'      ! Input PFT file logical name
      CHARACTER*16  :: LAITXTF = 'LAITXTF'      ! Input LAI file logical name
      CHARACTER*16  :: EFMAPS = 'EFMAPS'     ! Output EF file logical name
      CHARACTER*16  :: PFTS16 = 'PFTS16'     ! Output PFT file logical name
      CHARACTER*16  :: LAIS46 = 'LAIS46'     ! Output LAI file logical name

      CHARACTER*16  :: GDNAM

      CHARACTER*16  :: PROGNAME = 'TXT2IOAPI'

!***********************************************************************
!... Begin program
!... Initialize log file unit
      LDEV = INIT3() 

!... Get output file grid parameters
      MESG = 'Coordinate name: '
      CALL ENVSTR( 'GDNAM3D', MESG, 'ASACA36km', GDNAM, IOS )

!... Convert EFS
      RUN_EFS = ENVYN( 'RUN_EFS','Run TXT to IOAPI for EFS',
     &                .TRUE.,IOS)
      IF (IOS .GT. 0) THEN
        STOP 'Bad value for RUN_EFS'
      ENDIF
      IF (RUN_EFS) THEN
        CALL CNVT_EFS(GDNAM,EFSTXTF,EFMAPS)
      ENDIF

!... Convert PFT
      RUN_PFT = ENVYN( 'RUN_PFT','Run TXT to IOAPI for PFT',
     &                .TRUE.,IOS)
      IF (IOS .GT. 0) THEN
        STOP 'Bad value for RUN_PFT'
      ENDIF
      IF (RUN_PFT) THEN
        CALL CNVT_PFT(GDNAM,PFTTXTF,PFTS16)
      ENDIF

!... Convert LAI
      RUN_LAI = ENVYN( 'RUN_LAI','Run TXT to IOAPI for LAI',
     &                .TRUE.,IOS)
      IF (IOS .GT. 0) THEN
        STOP 'Bad value for RUN_LAI'
      ENDIF
      IF (RUN_LAI) THEN
        CALL CNVT_LAI(GDNAM,LAITXTF,LAIS46)
      ENDIF

!... Exit and close file
      CALL M3EXIT(PROGNAME,0,0,' ',0)

      ENDPROGRAM
