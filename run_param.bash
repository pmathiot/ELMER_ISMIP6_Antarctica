#!/bin/bash
# WARNING : this script only allows run witch begins in the 01/01 

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
WALLTIME=50000
NSTEP=365
TIME_STP=5 # in days

# first year in atmospheric forcing file / first year to read in the simulation 
START_YEAR_FORCING=1995    
START_SIMU=2015
OFFSET=$((START_SIMU-START_YEAR_FORCING))

# first year in oceanic forcing file / first year to read in the simulation
START_YEAR_FORCING_OC=1995
OFFSETOC=$((START_SIMU-START_YEAR_FORCING_OC))

calc() { awk "BEGIN{print $*}"; }
TIME_RST=`calc $NSTEP*$TIME_STP` # in days

# restart path and rst file (assume all in $IELMER)
RSTINITpath=${IELMER}/RST_simplified2/
RSTINITfile=ISMIP6_C-histAE_020.restart_indays1.nc

# MSH path and file
MSHINITpath=${IELMER}/MSH_simplified2/
