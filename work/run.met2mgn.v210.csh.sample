#!/bin/csh
#
# MET2MGN v2.10 
# --
#
#
# TPAR2IOAPI v2.03a 
# --added 26-category landuse capability for mm5camx (number of landuse categories defined by NLU) 
# --added capability for LATLON and UTM projections
# --added capability for MCIP v3.3 input (2m temperatures)
# --bug in PAR processing subroutine fixed where first few hours in GMT produced zero PAR
# --added code to fill missing par data (if valid data exists for the hours surrounding it)
#
# TPAR2IOAPI v2.0
# --added capability for MM5 or MCIP input
# 
#
#        RGRND/PAR options:
#           setenv MM5RAD  Y   Solar radiation obtained from MM5
#           OR 
#           setenv MCIPRAD Y   Solar radiation obtained from MCIP
#                  --MEGAN will internally calculate PAR for each of these options and user needs to  
#                    specify `setenv PAR_INPUT N' in the MEGAN runfile
#           OR
#           setenv SATPAR Y (satellite-derived PAR from UMD GCIP/SRB files)
#                  --user needs to specify `setenv PAR_INPUT Y' in the MEGAN runfile
#
#        TEMP options:
#           setenv CAMXTEMP Y         2m temperature, calculated from mm5camx output files
#           OR
#           setenv MM5MET  Y         2m temperature, calculated from MM5 output files
#                                     Note: 2m temperature is calculated since the P-X/ACM PBL
#                                     MM5 configuration (most commonly used LSM/PBL scheme for AQ 
#                                     modeling purposes) does not produce 2m temperatures.
#           OR
#           setenv MCIPMET Y         temperature obtained from MCIP
#              -setenv TMCIP  TEMP2   2m temperature, use for MCIP v3.3 or newer
#              -setenv TMCIP  TEMP1P5 1.5m temperature, use for MCIP v3.2 or older
#
#        TZONE   time zone for input mm5CAMx files 
#        NLAY    number of layers contained in input mm5CAMx files 
#        NLU     number of landuse categories contained in CAMx landuse file 
#

############################################################



############################################################
# Episodes
############################################################
set dom = 36 
set STJD = 2008201
set EDJD = 2008201

setenv EPISODE_SDATE 2008201
setenv EPISODE_STIME  000000    

############################################################
#set for grid
############################################################
setenv GRIDDESC GRIDDESC
setenv GDNAM3D USA${dom}km


############################################################
# Setting up directories and common environment variable
############################################################
source /data3/home/xjiang/MEGAN/MEGANv2.10/setcase.csh

setenv PROG met2mgn
setenv EXE $MGNEXE/$PROG


set logdir = logdir/$PROG
if ( ! -e $logdir) mkdir -p $logdir

set INPPATH     = /data3/home/xjiang/MCIP/mcip_out
set OUTPATH     = $MGNINP/MGNMET
if (! -e $OUTPATH) mkdir $OUTPATH


setenv PFILE $OUTPATH/PFILE
rm -fv $PFILE

############################################################
# Looping
############################################################
set JDATE = $STJD
set Y4 = 2008
set Y2 = 08 
set MM = 07
 DD = 21 
 DDm1 = 20
while ($JDATE <= $EDJD)


if ($JDATE == 2008367) set JDATE = 2009001
@ jdy  = $JDATE - 2000000
#set Y4 = `j2g $JDATE | awk '{print $1}'`
#set Y2 = `echo $Y4 | cut -c 3-4`
#set MM = `j2g $JDATE | awk '{print $2}'`
#set DD = `j2g $JDATE | awk '{print $3}'`

set Y4 = 2008
set Y2 = 08 
set MM = 07
@ DD++

@ JDATEm1 = $JDATE - 1
if ($JDATEm1 == 2008000) set JDATEm1 = 2007365
@ jdym1  = $JDATEm1 - 2000000
#set Y4m1 = `j2g $JDATEm1 | awk '{print $1}'`
#set Y2m1 = `echo $Y4m1 | cut -c 3-4`
#set MMm1 = `j2g $JDATEm1 | awk '{print $2}'`
#set DDm1 = `j2g $JDATEm1 | awk '{print $3}'`
set Y4m1 = 2008
set Y2m1 = 08
set MMm1 = 07
@ DDm1++

#set start/end dates
setenv STDATE ${jdy}00
setenv ENDATE ${jdy}23

#TEMP/PAR input choices
#
#set if using MM5 output files
setenv MM5MET N
setenv MM5RAD N
#setenv numMM5 2
#setenv MM5file1 /pete/pete5/fcorner/met/links/MMOUT_DOMAIN1_G$Y4$MM$DD
#setenv MM5file2 /pete/pete5/fcorner/met/links/MMOUT_DOMAIN1_G$Y4$MM$DD

#set if using UMD satellite PAR data
set PARDIR = $MGNINP/PAR
setenv SATPAR N
set satpar1 = "$PARDIR/$Y2m1${MMm1}par.h"
set satpar2 = "$PARDIR/$Y2${MM}par.h"

if ($satpar1 == $satpar2) then
  setenv numSATPAR 1
  setenv SATPARFILE1 $satpar2
else
  setenv numSATPAR 2
  setenv SATPARFILE1 $satpar1
  setenv SATPARFILE2 $satpar2
endif

#set if using MCIP output files
setenv MCIPMET Y
setenv TMCIP  TEMP2          #MCIP v3.3 or newer
#setenv TMCIP  TEMP1P5       #MCIP v3.2 or older

setenv MCIPRAD Y 
if ($JDATE == $EPISODE_SDATE) then
  setenv METCRO2Dfile1 $INPPATH/METCRO2D.$GDNAM3D
else
  setenv METCRO2Dfile1 $INPPATH/METCRO2D.$GDNAM3D
  setenv METCRO2Dfile2 $INPPATH/METCRO2D.$GDNAM3D
endif
setenv METCRO3Dfile  $INPPATH/METCRO3D.$GDNAM3D
setenv METDOT3Dfile  $INPPATH/METDOT3D.$GDNAM3D

setenv OUTFILE $OUTPATH/MET.MEGAN.$GDNAM3D.$JDATE.ncf
rm -rf $OUTFILE

$EXE |tee $logdir/log.$PROG.$GDNAM3D.$JDATE.txt 

@ JDATE++
end  # End while JDATE
