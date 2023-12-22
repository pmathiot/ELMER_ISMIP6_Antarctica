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
ENDITER=1
# define length of each segments
WALLTIME=12000
NSTEP=438
TIME_STP=2.5 # in days
calc() { awk "BEGIN{print $*}"; }
TIME_RST=`calc $NSTEP*$TIME_STP` # in days

#------------------------------------------------------------------------------
#                               FORCING DATA
#------------------------------------------------------------------------------
# first year to read in the simulation 
START_SIMU=2015

# time origine for output
TIME_ORIGINE=1995

# nature of oceanic forcing
# constant : for constant forcing
# variable : for variable forcing
FORCING='constant'

# smb without or with anomaly (constant, anomaly or variable)
# constant : forcing is read at the beginning of the simulation only
# anomaly : sum of a constant forcing and an anomaly forcing which varie each year
# variable : forcing read each year
SMB_METHOD='constant'

# The starting year of forcing files is detected automatically by
# prepare_elmer.bash. To override this automatic detection, you can set one or
# both of the following variables (for the atmospheric and oceanic forcing
# files, respectively)
# first year in atmospheric forcing file
# START_YEAR_FORCING=1995
# first year in oceanic forcing file
# START_YEAR_FORCING_OC=1995

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
# friction (linear, weertman or regularized coulomb)
FRICTION='linear'

# melt parameterisation (PICO or QUADRATIC)
# for QUADRATIC : precise local or global for the calculation of slope
# for QUADRATIC : precise correction : TRUE or FALSE for adding T correction of oceanic temperature
MELT='QUADRATIC'
SLOPE='local'
CORRECTION='FALSE'

