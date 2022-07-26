#!/bin/bash
######################
## IRENE   TGCC/CEA ##
######################
#MSUB -r interp
#MSUB -o interp.o%j
#MSUB -e interp.e%j
#MSUB -n 1 
#MSUB -x (eq -exclusive, remove it for share queue)
#MSUB -T 90 (walltime in second)
#MSUB -A gen6066 (gen6035 for us)
#MSUB -q rome (partition name, rome for us)
#MSUB -m work,scratch (spaces accessible from the compute node)

set -x

module load cdo
newgrp  gen6066

cd ${BRIDGE_MSUB_PWD}

export OMP_NUM_THREADS=32

## conservative weights from ISM_ANT_4km to elmer ant50.gl1 grid
cdo genycon,ant50.gl1-ismip6_grid.nc racmo_wessem_initmip8km_mean1995_2014_smbice.nc yCONweights_racmo_ant50.gl1.nc
