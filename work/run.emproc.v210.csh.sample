#! /bin/csh -f
########################################################################
source /data3/home/xjiang/MEGAN/MEGANv2.10/setcase.csh
## Directory setups
setenv PRJ USA 
setenv PROMPTFLAG N

# Program directory
setenv PROG   emproc
setenv EXEDIR $MGNEXE
setenv EXE    $EXEDIR/$PROG

# Input map data directory
setenv INPDIR $MGNINP/MAP

# MCIP input directory
setenv METDIR $MGNINP/MGNMET

# Intermediate file directory
setenv INTDIR $MGNINT

# Output directory
setenv OUTDIR $MGNOUT

# Log directory
setenv LOGDIR $MGNLOG/$PROG
if ( ! -e $LOGDIR ) mkdir -p $LOGDIR
########################################################################

set dom = 36 
set JD = 2008201
while ($JD <= 2008201)
########################################################################
# Set up time and date to process
setenv SDATE $JD        #start date
setenv STIME 0
setenv RLENG 240000
########################################################################

########################################################################
# Set up for MEGAN
setenv RUN_MEGAN   Y       # Run megan?


# By default MEGAN will use data from MGNMET unless specify below
setenv ONLN_DT     Y       # Use online daily average temperature
                           # No will use from EFMAPS

setenv ONLN_DS     Y       # Use online daily average solar radiation
                           # No will use from EFMAPS

# Grid definition
setenv GRIDDESC $MGNRUN/GRIDDESC
setenv GDNAM3D USA${dom}km 

# EFMAPS
setenv EFMAPS $INPDIR/EFMAPS.${PRJ}${dom}.ncf

# PFTS16
setenv PFTS16 $INPDIR/PFTS16.${PRJ}${dom}.ncf

# LAIS46
setenv LAIS46 $INPDIR/LAIS46.${PRJ}${dom}.ncf

# MGNMET
setenv MGNMET $METDIR/MET.MEGAN.$GDNAM3D.${SDATE}.ncf

# Output
setenv MGNERS $INTDIR/ER.$GDNAM3D.${SDATE}.ncf

########################################################################
## Run MEGAN
if ( $RUN_MEGAN == 'Y' ) then
   rm -f $MGNERS
   time $EXE >&! $LOGDIR/log.run.$PROG.$GDNAM3D.$SDATE.txt
endif

@ JD++
end  # End while JD
