#!/bin/bash
# WARNING : this script only allows run witch begins in the 01/01 

#----------------------------------------------------------------------------
#                                 COMPUTING DATA
#----------------------------------------------------------------------------
# number of elmer partition + xios server (48 + 1 in this case)
NELMER=48
NXIOS=1
NP=$((NELMER+NXIOS))

# number of HPC nodes
NN=1

# first iteration number (if more than 1, means restart i-1 are already in place)
# end iteration
STARTITER=1
ENDITER=2

# define length of each segments
WALLTIME=12000
NSTEP=146
TIME_STP=2.5 # in days
calc() { awk "BEGIN{print $*}"; }
TIME_RST=`calc $NSTEP*$TIME_STP` # in days

#------------------------------------------------------------------------------
#                               FORCING DATA
#------------------------------------------------------------------------------
# first year in atmospheric forcing file / first year to read in the simulation 
START_YEAR_FORCING=1995    
START_SIMU=1995
OFFSET=$((START_SIMU-START_YEAR_FORCING))

# first year in oceanic forcing file
START_YEAR_FORCING_OC=1995
OFFSETOC=$((START_SIMU-START_YEAR_FORCING_OC))

# time orginie for output
TIME_ORIGINE=1995
#------------------------------------------------------------------------------
#                               RESTART DATA
#------------------------------------------------------------------------------
# restart path and rst file (assume all in $IELMER)
RSTINITpath=${IELMER}/RST_simplified2/
RSTINITfile=ISMIP6_W1-HistAE_002.restart_frc.nc 

#------------------------------------------------------------------------------
#                               MESH DATA
#------------------------------------------------------------------------------
# MSH path and file
MSHINITpath=${IELMER}/MSH_simplified2/

#------------------------------------------------------------------------------
#                               PARAM DATA
#------------------------------------------------------------------------------
# smb without or with anomaly (constant, anomaly or variable)
# constant : forcing is read at the beginning of the simulation only
# anomaly : sum of a constant forcing and an anomaly forcing which varie each year
# variable : forcing read each year
SMB_METHOD='anomaly'

# friction (linear, weertman or regularized coulomb)
FRICTION='linear'

# melt parameterisation (PICO or QUADRATIC)
# for QUADRATIC : precise local or global for the calculation of slope
# for QUADRATIC : precise correction : TRUE or FALSE for adding T correction of oceanic temperature
MELT='QUADRATIC'
SLOPE='local'
CORRECTION='FALSE'

