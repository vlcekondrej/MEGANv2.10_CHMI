#! /bin/csh -f
########################################################################
source /data3/home/xjiang/MEGAN/MEGANv2.10/setcase.csh
## Directory setups
setenv PRJ USA 
setenv PROMPTFLAG N

# Program directory
setenv PROG   mgn2mech.wmap
setenv EXEDIR $MGNEXE
setenv EXE    $EXEDIR/$PROG

# Input map data directory
setenv INPDIR $MGNINP/MAP

# Intermediate file directory
setenv INTDIR $MGNOUT/INT

# Output directory
setenv OUTDIR $MGNOUT

# MCIP input directory
setenv METDIR $MGNINP/MGNMET

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
setenv TSTEP 10000
########################################################################

########################################################################
# Set up for MECHCONV
setenv RUN_SPECIATE   Y    # run MG2MECH

setenv RUN_CONVERSION Y    # run conversions?
                           # run conversions MEGAN to model mechanism
                           # units are mole/s

setenv SPCTONHR       N    # speciation output unit in tonnes per hour
                           # This will convert 138 species to tonne per
                           # hour or mechasnim species to tonne per hour.
                           
# If RUN_CONVERSION is set to "Y", one of mechanisms has to be selected.
setenv MECHANISM    RADM2
#setenv MECHANISM    RACM
#setenv MECHANISM    CBMZ
#setenv MECHANISM    CB05
#setenv MECHANISM    CB6
#setenv MECHANISM    SOAX
#setenv MECHANISM    SAPRC99
#setenv MECHANISM    SAPRC99Q
#setenv MECHANISM    SAPRC99X

# Grid name
setenv GDNAM3D USA${dom}km 

# EFMAPS NetCDF input file
setenv EFMAPS  $INPDIR/EFMAPS.${PRJ}${dom}.ncf

# PFTS16 NetCDF input file
setenv PFTS16  $INPDIR/PFTS16.${PRJ}${dom}.ncf

# MEGAN ER filename
setenv MGNERS $INTDIR/ER.$GDNAM3D.${SDATE}.ncf

# Output filename
setenv MGNOUT $OUTDIR/MEGANv2.10.$GDNAM3D.$MECHANISM.$SDATE.ncf

########################################################################
## Run speciation and mechanism conversion
if ( $RUN_SPECIATE == 'Y' ) then
   rm -f $MGNOUT
   $EXE | tee $LOGDIR/log.run.$PROG.$GDNAM3D.$MECHANISM.$SDATE.txt
endif

@ JD++
end  # End while JD
