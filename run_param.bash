#!/bin/bash

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
NSTEP=10
TIME_STP=1.0       # in days

calc() { awk "BEGIN{print $*}"; }
TIME_RST=`calc $NSTEP*$TIME_STP` # in days

# restart path and rst file (assume all in $IELMER)
RSTINITpath=${IELMER}/RST_simplified2/
RSTINITfile=restart_newmesh_beta_coulomg_reg_time_eq_0.nc

# MSH path and file
MSHINITpath=${IELMER}/MSH_simplified2/
