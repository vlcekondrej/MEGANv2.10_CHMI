
      SUBROUTINE CNVT_EFS(GDNAM,EFSTXTF,EFMAPS)
!***********************************************************************
!  This subroutine converts EFS txt to EFS IOAPI
!
!  Called by TXT2IOAPI
!
!  Call:
!     PARSLINE
!     CHECKMEM
!
!  Created by Tan  07/07/2011 for MEGAN v2.10
!
!  History:
!***********************************************************************

      IMPLICIT NONE

C...........   INCLUDES:
      INCLUDE 'PARMS3.EXT'    !  I/O API parameters
      INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
      INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures

      INCLUDE 'SPC_MGN.EXT'   !  MEGAN species

C...........   EXTERNAL FUNCTIONS and their descriptions:
      INTEGER      PROMPTFFILE, STR2INT
      REAL         STR2REAL
      LOGICAL      DSCGRID
      INTEGER      INDEX1
      EXTERNAL     PROMPTFFILE, DSCGRID
     &             STR2INT, STR2REAL, INDEX1

C.........  LOCAL VARIABLES and their descriptions:
      INTEGER       IOS            ! i/o status
      INTEGER       ILINE          ! current line
      INTEGER       ADEV           ! unit no. for input data file
      INTEGER       I, J, S, N     ! Counters
      CHARACTER*16    CNAME        ! Coord name
      CHARACTER(LEN=10000) LINE     ! input ine buffer
      CHARACTER*256    MESG        ! message buffer
      INTEGER       CID, INX, INY  ! Input grid x and y

      CHARACTER*16  :: EFSTXTF
      CHARACTER*16  :: EFMAPS
      CHARACTER*16  :: GDNAM

      ! CID, I, J, LAT, LONG, SR, TEMP, 20 EFS
      INTEGER, PARAMETER :: MXTCOL = 27      ! Columns in an input line
      CHARACTER*30     SEGMENT( MXTCOL )     ! Input line fields
      INTEGER          XREF(MXTCOL)

      INTEGER, PARAMETER :: NVARS = 24       ! Number of output variables
      CHARACTER*16 VNAME( NVARS ) ! Output variables

      REAL, ALLOCATABLE :: IBUFF( :, :, : )  ! input, output buffer

      CHARACTER*16  :: PROGNAME = 'CNVT_EFS'

!***********************************************************************
!... Begin program

!... Get output file grid parameters
      GDNAM3D = GDNAM
      IF( .NOT. DSCGRID( GDNAM3D, CNAME, GDTYP3D,
     &              P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,
     &              XORIG3D, YORIG3D, XCELL3D, YCELL3D,
     &              NCOLS3D, NROWS3D, NTHIK3D ) ) THEN
         MESG = 'Could not get grid description.'
         CALL M3EXIT ( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      print*,'NCOLS3D ',GDNAM3D, CNAME, GDTYP3D,
     &              P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,
     &              XORIG3D, YORIG3D, XCELL3D, YCELL3D,
     &              NCOLS3D, NROWS3D, NTHIK3D

!... Allocate memory
      ALLOCATE ( IBUFF( NVARS, NCOLS3D, NROWS3D ),
     &           STAT = IOS )
      CALL CHECKMEM ( IOS, 'IBUFF', PROGNAME )
      IBUFF = 1.0

!... Get input file unit
      MESG = 'Enter logical name of EFS data file.'
      ADEV = PROMPTFFILE( MESG, .TRUE., .TRUE., EFSTXTF, PROGNAME )
      IF( ADEV < 0 ) THEN
         MESG = 'Error opening input EFS file.'
         CALL M3EXIT ( PROGNAME, 0, 0, MESG, 2 )
      ENDIF

!... Set VNAME
      DO N = 1,N_MGN_SPC
        VNAME( N+4 ) = 'EF_'//TRIM(MGN_SPC(N))
      ENDDO

!... Read TXT input file
      ILINE = 0
      DO
         READ( ADEV, 1000, IOSTAT = IOS ) LINE
         ILINE = ILINE + 1
         IF( IOS < 0 ) EXIT       ! End of file, exit loop
         IF( LINE == ' ' ) CYCLE  ! Skip blank line
         IF( IOS > 0 ) THEN
            WRITE( MESG, 1010 ) 'I/O error', IOS,
     &            'reading input file at line', ILINE
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
         ENDIF

         CALL PARSLINE( LINE, MXTCOL, SEGMENT )  ! break lines into segments
         IF ( ILINE .EQ. 1 ) THEN
            ! Assumption
            ! Column 1 : CID
            ! Column 2 : ICELL
            ! Column 3 : JCELL
            ! Column 4 : LAT
            ! Column 5 : LONG
            ! column 6 : SR
            ! Column 7 : TEMP
            DO S = 1 , MXTCOL - 3
               MESG = SEGMENT( S + 3 )
               CALL UPCASE( MESG )  ! to all caps
               IF (       TRIM( MESG ) .EQ. 'LAT'    ) THEN
                  XREF( S ) = 1
                  VNAME( XREF(S) ) = TRIM( MESG )
               ELSE IF (  TRIM( MESG ) .EQ. 'LONG'   ) THEN
                  XREF( S ) = 2
                  VNAME( XREF(S) ) = TRIM( MESG )
               ELSE IF (  TRIM( MESG ) .EQ. 'DSRAD'  ) THEN
                  XREF( S ) = 3
                  VNAME( XREF(S) ) = TRIM( MESG )
               ELSE IF (  TRIM( MESG ) .EQ. 'DTEMP'  ) THEN
                  XREF( S ) = 4
                  VNAME( XREF(S) ) = TRIM( MESG )
               ELSE IF ( INDEX1(TRIM( MESG ),N_MGN_SPC ,MGN_SPC)
     &                                                      .NE. 0) THEN
                  N = INDEX1(TRIM( MESG ),N_MGN_SPC ,MGN_SPC)
                  XREF( S ) = N + 4
                  VNAME( XREF(S) ) = 'EF_'//TRIM( MESG )
               ELSE
                  write(*,*)MESG//' is not valid'
                  STOP
               ENDIF

               PRINT*,'READ VNAME:',VNAME(XREF(S))
            ENDDO
            CYCLE  ! go back to read 2nd line
         ENDIF

         CID = STR2INT( SEGMENT(1) )   ! convert character to integer
         INX = STR2INT( SEGMENT(2) )   ! convert character to integer
         INY = STR2INT( SEGMENT(3) )   ! convert character to integer
        print*,'CID ',cid,inx,iny,xref
         DO S = 1, 4
            IBUFF(S,INX,INY) = STR2REAL( SEGMENT( 3+S ) )        ! convert char to real
         print*,'IBUFF ',S,IBUFF(S,INX,INY)
         ENDDO
         DO S = 5, MXTCOL - 3
            IBUFF(XREF(S),INX,INY) = STR2REAL( SEGMENT( 3+S ) )  ! convert char to real
         print*,'IBUFF ',S,IBUFF(S,INX,INY) 
         ENDDO
      ENDDO

      WRITE( MESG, 1010 ) 'Total number of input lines read:',ILINE
      CALL M3MESG( MESG )

!... Setup additional outfile parameters and open outfile
      FTYPE3D = GRDDED3
      NVARS3D = NVARS
      VNAME3D(1:NVARS) = VNAME(1:NVARS)
      VTYPE3D(1:NVARS) = M3REAL
      TSTEP3D = 0 ! time independent
      NLAYS3D = 1
      CALL NAMEVAL ( EFSTXTF, MESG )  ! get input file name
      FDESC3D( 1 ) = 'Converted from MEGAN input text file: '
      FDESC3D( 2 ) = TRIM(MESG)
      DO S = 1, NVARS
! For LAT, and LONG
         IF ( INDEX( VNAME3D( S ), 'LAT'      ) .EQ. 1 )
     &                                 UNITS3D( S ) = 'DEGREE'
         IF ( INDEX( VNAME3D( S ), 'LONG'     ) .EQ. 1 )
     &                                 UNITS3D( S ) = 'DEGREE'
! For daily SRAD, and daily TEMP
         IF ( INDEX( VNAME3D( S ), 'D_SRAD'   ) .EQ. 1 )
     &                                 UNITS3D( S ) = 'Watt/m2'
         IF ( INDEX( VNAME3D( S ), 'D_TEMP'   ) .EQ. 1 )
     &                                 UNITS3D( S ) = 'K'
! For MEGAN species
         DO N = 1, N_MGN_SPC
            IF ( INDEX( VNAME3D( S ), 'EF_'//MGN_SPC( N ) ) .EQ. 1 )
     &                     UNITS3D( S ) = 'ug'//MGN_SPC( N )//'/m2.hr'
         ENDDO

      ENDDO
      IF ( .NOT. OPEN3( EFMAPS, FSCREA3, PROGNAME ) ) THEN
           MESG = 'Could not open file "'//TRIM(EFMAPS)//'" for output'
           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF

!... Write output
      DO S = 1, NVARS
         PRINT*,'min cell value for:',TRIM(VNAME3D(S)),MINVAL(IBUFF(S,:,:))
         PRINT*,'max cell value for:',TRIM(VNAME3D(S)),MAXVAL(IBUFF(S,:,:))
         IF ( .NOT. WRITE3(EFMAPS,VNAME3D(S),0,0,IBUFF(S,:,:))) THEN
            MESG = 'Error writing to  file "'//TRIM(EFMAPS)//'"'
            CALL M3EXIT(PROGNAME,0,0,MESG,2)
         ENDIF
      ENDDO

      IF ( CLOSE3(EFMAPS) ) THEN
        CALL M3MESG( 'Closing EFS IOAPI file')
      ENDIF
      CLOSE(ADEV)

1000  FORMAT( 40A )
1010  FORMAT( 40( A, :, I8, :, 1X ) )

      END SUBROUTINE
