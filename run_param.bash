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
ENDITER=5

# define length of each segments
WALLTIME=10800
NSTEP=365
TIME_STP=1 # in days

# first year in forcing file / first year to read in the simulation 
START_YEAR_FORCING=1995    
START_SIMU=1995
OFFSET=$((START_SIMU-START_YEAR_FORCING))

calc() { awk "BEGIN{print $*}"; }
TIME_RST=`calc $NSTEP*$TIME_STP` # in days

# restart path and rst file (assume all in $IELMER)
RSTINITpath=${IELMER}/RST_simplified2/
RSTINITfile=restart_newmesh_t0_beta_coulomg_reg.nc

# MSH path and file
MSHINITpath=${IELMER}/MSH_simplified2/
