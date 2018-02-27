      LOGICAL FUNCTION CKGDIOAPI2D( INFILE1, INFILE2 )

!***********************************************************************
!
!  DESCRIPTION:
!      This function compares two ioapi file for
!       - grid/projection parameters
!       - TRUE if all match
!
!  PRECONDITIONS REQUIRED:
!       - INFILE1 and INFILE2 are already opened
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!    started 10/04 by Jack Chen prototype
!
!**************************************************************************
!...........   MODULES for public variables
      IMPLICIT NONE

!...........   INCLUDES:
      INCLUDE 'PARMS3.EXT'    !  I/O API parameters
      INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
      INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

!............  External Functions:
!     CHARACTER*16   PROMPTMFILE
!     LOGICAL        CHKMETEM
!     EXTERNAL       PROMPTMFILE, CHKMETEM
              
!........... Argument variables:
      CHARACTER(LEN=*), INTENT (IN) :: INFILE1 
      CHARACTER(LEN=*), INTENT (IN) :: INFILE2

!........... Local Variables
      INTEGER :: NCOLS1 = 0      ! number of columns in grid
      INTEGER :: NGRID1 = 0      ! number of cells in grid
      INTEGER :: NLAYS1 = 0      ! number of layers
      INTEGER :: NROWS1 = 0      ! number of rows in grid
      INTEGER :: NTHIK1 = 1
      INTEGER :: GDTYP1 = -1     ! i/o api grid type code
      REAL    :: P_ALP1 = 0.D0   ! projection alpha
      REAL    :: P_BET1 = 0.D0   ! projection beta
      REAL    :: P_GAM1 = 0.D0   ! projection gamma
      REAL    :: XCENT1 = 0.D0   ! x-center of projection
      REAL    :: YCENT1 = 0.D0   ! y-center of projection
      REAL    :: XORIG1 = 0.D0   ! x-origin of grid
      REAL    :: YORIG1 = 0.D0   ! y-origin of grid
      REAL    :: XCELL1 = 0.D0   ! x-dim of cells
      REAL    :: YCELL1 = 0.D0   ! y-dim of cells
      INTEGER :: VGTYP1 = -1     ! type of vertical coordinates
      REAL    :: VGTOP1 = 0.0    ! model-top, for sigma coord types
      REAL,ALLOCATABLE,DIMENSION(:)  :: VGLVS1     ! vertical coordinate values

      INTEGER :: NCOLS2 = 0      ! number of columns in grid
      INTEGER :: NGRID2 = 0      ! number of cells in grid
      INTEGER :: NLAYS2 = 0      ! number of layers
      INTEGER :: NROWS2 = 0      ! number of rows in grid
      INTEGER :: NTHIK2 = 1
      INTEGER :: GDTYP2 = -1     ! i/o api grid type code
      REAL    :: P_ALP2 = 0.D0   ! projection alpha
      REAL    :: P_BET2 = 0.D0   ! projection beta
      REAL    :: P_GAM2 = 0.D0   ! projection gamma
      REAL    :: XCENT2 = 0.D0   ! x-center of projection
      REAL    :: YCENT2 = 0.D0   ! y-center of projection
      REAL    :: XORIG2 = 0.D0   ! x-origin of grid
      REAL    :: YORIG2 = 0.D0   ! y-origin of grid
      REAL    :: XCELL2 = 0.D0   ! x-dim of cells
      REAL    :: YCELL2 = 0.D0   ! y-dim of cells
      INTEGER :: VGTYP2 = -1     ! type of vertical coordinates
      REAL    :: VGTOP2 = 0.0    ! model-top, for sigma coord types
      REAL,ALLOCATABLE,DIMENSION(:) :: VGLVS2     ! vertical coordinate values

      INTEGER STATUS
      INTEGER L
      CHARACTER*256   MESG
      CHARACTER*16 :: PROGNAME = 'CKGDIOAPI2D' ! Program name

!........... Begin subroutine

!... Open 1st input file and store variables
      IF( .NOT. DESC3( INFILE1 ) ) THEN
      MESG = 'Could not get input file "'//TRIM(INFILE1) 
     &    //'" description'
          CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF

      GDTYP1 = GDTYP3D
      P_ALP1 = P_ALP3D
      P_BET1 = P_BET3D
      P_GAM1 = P_GAM3D
      XCENT1 = XCENT3D
      YCENT1 = YCENT3D
      XORIG1 = XORIG3D
      YORIG1 = YORIG3D
      NROWS1 = NROWS3D
      NCOLS1 = NCOLS3D
      XCELL1 = XCELL3D
      YCELL1 = YCELL3D
      NTHIK1 = NTHIK3D
      NLAYS1 = NLAYS3D
      IF( NLAYS1 .GT. 1 ) THEN
         ALLOCATE ( VGLVS1 ( NLAYS1 + 1 ), STAT = STATUS )
         CALL CHECKMEM( STATUS, 'VGLVS1', PROGNAME )
         VGTYP1 = VGTYP3D
         VGLVS1 = VGLVS3D  ! array
         VGTOP1 = VGTOP3D
      ENDIF

!... Open 2nd input file and store variables
      IF( .NOT. DESC3( INFILE2 ) ) THEN
          MESG = 'Could not get input file "'//TRIM(INFILE2)// 
     &    '" description' 
          CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF

      GDTYP2 = GDTYP3D
      P_ALP2 = P_ALP3D
      P_BET2 = P_BET3D
      P_GAM2 = P_GAM3D
      XCENT2 = XCENT3D
      YCENT2 = YCENT3D
      XORIG2 = XORIG3D
      YORIG2 = YORIG3D
      NROWS2 = NROWS3D
      NCOLS2 = NCOLS3D
      XCELL2 = XCELL3D
      YCELL2 = YCELL3D
      NTHIK2 = NTHIK3D
      NLAYS2 = NLAYS3D
      IF( NLAYS2 .GT. 1 ) THEN
         ALLOCATE ( VGLVS2 ( NLAYS2 + 1 ), STAT = STATUS )
         CALL CHECKMEM( STATUS, 'VGLVS2', PROGNAME )
         VGTYP2 = VGTYP3D
         VGLVS2 = VGLVS3D  ! array
         VGTOP2 = VGTOP3D
      ENDIF

      IF( GDTYP1 .NE. GDTYP2 ) THEN
          WRITE( MESG, * )'Differences in GDTYP: ',GDTYP1,' and ',GDTYP2
          CALL M3WARN( PROGNAME, 0, 0, MESG ) 
          CKGDIOAPI2D = .FALSE.
          RETURN
      ENDIF

      IF( NROWS1 .NE. NROWS2 ) THEN
          WRITE( MESG, * )'Differences in NROWS: ',NROWS1,' and ',NROWS2
          CALL M3WARN( PROGNAME, 0, 0, MESG ) 
          CKGDIOAPI2D = .FALSE.
          RETURN
      ENDIF

      IF( NCOLS1 .NE. NCOLS2 ) THEN 
          WRITE( MESG, * )'Differences in NCOLS: ',NCOLS1,' and ',NCOLS2
          CALL M3WARN( PROGNAME, 0, 0, MESG ) 
          CKGDIOAPI2D = .FALSE.
          RETURN
      ENDIF

      IF( FLTERR ( P_ALP1 , P_ALP2 ) ) THEN
          WRITE( MESG, * )'Differences in P_ALP: ',P_ALP1,' and ',P_ALP2
          CALL M3WARN( PROGNAME, 0, 0, MESG ) 
          CKGDIOAPI2D = .FALSE.
          RETURN
      ENDIF

      IF( FLTERR ( P_BET1 , P_BET2 ) ) THEN
          WRITE( MESG, * )'Differences in P_BET: ',P_BET1,' and ',P_BET2
          CALL M3WARN( PROGNAME, 0, 0, MESG ) 
          CKGDIOAPI2D = .FALSE.
          RETURN
      ENDIF

      IF( FLTERR ( P_GAM1 , P_GAM2 ) ) THEN
          WRITE( MESG, * )'Differences in P_GAM: ',P_GAM1,' and ',P_GAM2
          CALL M3WARN( PROGNAME, 0, 0, MESG ) 
          CKGDIOAPI2D = .FALSE.
          RETURN
      ENDIF

      IF( FLTERR ( XCENT1 , XCENT2 ) ) THEN
          WRITE( MESG, * )'Differences in XCENT: ',XCENT1,' and ',XCENT2
          CALL M3WARN( PROGNAME, 0, 0, MESG ) 
          CKGDIOAPI2D = .FALSE.
          RETURN
      ENDIF

      IF( FLTERR ( YCENT1 , YCENT2 ) ) THEN
          WRITE( MESG, * )'Differences in YCENT: ',YCENT1,' and ',YCENT2
          CALL M3WARN( PROGNAME, 0, 0, MESG ) 
          CKGDIOAPI2D = .FALSE.
          RETURN
      ENDIF

      IF( FLTERR ( XORIG1 , XORIG2 ) ) THEN
          WRITE( MESG, * )'Differences in XORIG: ',XORIG1,' and ',XORIG2
          CALL M3WARN( PROGNAME, 0, 0, MESG ) 
          CKGDIOAPI2D = .FALSE.
          RETURN
      ENDIF

      IF( FLTERR ( YORIG1 , YORIG2 ) ) THEN
          WRITE( MESG, * )'Differences in YORIG: ',YORIG1,' and ',YORIG2
          CALL M3WARN( PROGNAME, 0, 0, MESG ) 
          CKGDIOAPI2D = .FALSE.
          RETURN
      ENDIF

      IF( FLTERR ( XCELL1 , XCELL2 ) ) THEN
       WRITE( MESG, * )'Differences in XCELLS: ',XCELL1,' and ',XCELL2
          CALL M3WARN( PROGNAME, 0, 0, MESG ) 
          CKGDIOAPI2D = .FALSE.
          RETURN
      ENDIF

      IF( FLTERR ( YCELL1 , YCELL2 ) ) THEN
          WRITE( MESG, * )'Differences in YCELL: ',YCELL1,' and ',YCELL2
          CALL M3WARN( PROGNAME, 0, 0, MESG ) 
          CKGDIOAPI2D = .FALSE.
          RETURN
      ENDIF
      
      IF( NTHIK1 .NE. NTHIK2 ) THEN
          WRITE( MESG, * )'Differences in NTHIK: ',NTHIK1,' and ',NTHIK2
          CALL M3WARN( PROGNAME, 0, 0, MESG ) 
          CKGDIOAPI2D = .FALSE.
          RETURN
      ENDIF

      IF( NLAYS1 .NE. NLAYS2 ) THEN
          WRITE( MESG, * )'Differences in NLAYS: ',NLAYS1,' and ',NLAYS2
          CALL M3WARN( PROGNAME, 0, 0, MESG ) 
          CKGDIOAPI2D = .TRUE.
          RETURN
      ENDIF

      IF( NLAYS1 .GT. 1 ) THEN
         IF( VGTYP1 .NE. VGTYP2 ) THEN
       WRITE( MESG, * )'Differences in VGTYP: ',VGTYP1,' and ',VGTYP2
            CALL M3WARN( PROGNAME, 0, 0, MESG ) 
            CKGDIOAPI2D = .TRUE.
            RETURN
         ENDIF
          IF( VGTOP1  .NE. VGTOP2 ) THEN
       WRITE( MESG, * )'Differences in VGTOP: ',VGTOP1,' and ',VGTOP2
            CALL M3WARN( PROGNAME, 0, 0, MESG ) 
            CKGDIOAPI2D = .TRUE.
            RETURN
         ENDIF
         DO L = 1, NLAYS1+1
            IF( FLTERR ( VGLVS1( L ) , VGLVS2( L ) ) ) THEN
            WRITE( MESG, * )'Differences in VGLVS: ',VGLVS1( L ), 
     &  ' and ',VGLVS2( L )
               CALL M3WARN( PROGNAME, 0, 0, MESG ) 
               CKGDIOAPI2D = .TRUE.
               RETURN
            ENDIF
         ENDDO
      ENDIF ! NLAY

      CKGDIOAPI2D = .TRUE.
      
      RETURN

      CONTAINS

!... Internal Subprogram - comparing two real values
!... Copied from chkmetem.f of SMOKE program

      LOGICAL FUNCTION FLTERR( P, Q )

      REAL, INTENT (IN) ::  P, Q

      FLTERR =  
     &     ( (P - Q)**2  .GT.  1.0E-12*( P*P + Q*Q + 1.0E-5 ) )

      RETURN

      END FUNCTION FLTERR

      END FUNCTION CKGDIOAPI2D

