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
NSTEP=365
TIME_STP=5 # in days
calc() { awk "BEGIN{print $*}"; }
TIME_RST=`calc $NSTEP*$TIME_STP` # in days

#------------------------------------------------------------------------------
#                               FORCING DATA
#------------------------------------------------------------------------------
START_SIMU=2015
# The starting year of forcing files is detected automatically by
# prepare_elmer.bash. To override this automatic detection, you can set one or
# both of the following variables (for the atmospheric and oceanic forcing
# files, respectively)
# START_YEAR_FORCING=1995
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
MELT='PICO'
