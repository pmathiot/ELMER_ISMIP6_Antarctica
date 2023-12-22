#!/bin/bash
CONFIG=ANT50.GL1
CASE=ISMIP6

HELMER=`pwd`

# check if error fix or not
if [ -f zELMER_*ERROR* ]; then echo ' ERROR files still there. rm zERROR_*_* files if you fixed it before re-run'; exit 42; fi

# source arch file
. param_arch.bash

# load modules
load_elmer_modules

if [ ! -d WELMER  ]; then mkdir -p  $WELMER   ; fi
if [ ! -d RELMER  ]; then mkdir -p  $RELMER   ; fi
if [ ! -d SELMER  ]; then mkdir -p  $SELMER   ; fi

# source run param file
. run_param.bash


# to fix issue for vtu group gidbit
chgrp -R ${GROUPUSR} ${WELMER}

# set symbolic link
if [ ! -L MY_WORK    ]; then ln -s $WELMER MY_WORK    ; fi
if [ ! -L MY_RESTART ]; then ln -s $RELMER MY_RESTART ; fi
if [ ! -L MY_OUTPUT  ]; then ln -s $SELMER MY_OUTPUT  ; fi

echo ''
echo "Copy executable to $WELMER"
echo "=========================="
echo ''
# to clean this we could find a way to check only what is needed in the sif file
# like this workdir only contains what is needed
\cp MY_BLD/* $WELMER/.

# copy xios
\cp $XIOS_BIN/xios_server.exe $WELMER/.

# build app.conf file
touch app.conf
echo "$NELMER ElmerSolver_mpi" >> app.conf
echo "$NXIOS ./xios_server.exe" >> app.conf
mv app.conf $WELMER/.

# copy sif and param
echo ''
echo "Copy incf and param to: $WELMER"
echo "======================="
echo
\cp $CONFIG-${CASE}_elmer.param  $WELMER/elmer.param
\cp $CONFIG-${CASE}_elmer.incf   $WELMER/elmer.incf
\cp iodef_elmer.xml              $WELMER/.
#\cp context_elmer.xml            $WELMER/.
\cp iodef_pp.xml                 $WELMER/.
\cp context_pp.xml               $WELMER/.
\cp $CONFIG-${CASE}_pp.sif               $WELMER/pp.sif

# prepare file_def.xml
sed -e "s/<TIME_RST>/$TIME_RST/g"   ${CONFIG}-${CASE}_file_def.xml > $WELMER/file_def_elmer.xml

# link data
# elmer point directly to the data dir => no need to do anything
# except for the mesh
if [ ! -d $WELMER/MSH/partitioning.$NELMER ] ; then mkdir -p $WELMER/MSH/partitioning.$NELMER ; chmod g+s $WELMER/MSH ; chmod g+s $WELMER/MSH/partitioning.$NELMER ; fi
echo '   - link partionned mesh'
ln -sf $MSHINITpath/partitioning.$NELMER/part.* $WELMER/MSH/partitioning.$NELMER/.
echo '   - link serial mesh'
ln -sf $MSHINITpath/mesh.* $WELMER/MSH/.

echo ''
echo "$(($ENDITER - $STARTITER + 1)) segment to run"
echo '===================='

for ((ijob=$STARTITER ; ijob<=$ENDITER ; ijob++))
do
    im1=`printf %03d $((ijob-1))` ;
    ip1=`printf %03d $((ijob+1))` ;
    i=`printf %03d $((ijob))`     ;

    if [ -f zELMER_${i}_SUCCESSFUL ]; then echo "Run already successful. rm zELMER_${i}_SUCCESSFUL files if you want to overwrite this segment"; exit 42; fi
    
    # name
    NAME=$CONFIG-$CASE

    # get restart
    # no need to do anything, elmer point directly to the directory
    # except for the first one
    if [[ $((ijob-1)) -eq 0 ]]; then
      if [[ $RSTINITfile != NONE ]]; then
         RSTFILEnc=$RSTINITfile
         cp $RSTINITpath/${RSTINITfile} $WELMER/MSH/restart_${im1}.nc
         ncap2 -O -s "elmer_time=elmer_time*0" $WELMER/MSH/restart_${im1}.nc $WELMER/MSH/restart_${im1}.nc  
         #ln -sf $RSTINITpath/${RSTINITfile} $WELMER/MSH/restart_${im1}.nc
      fi
    else
         RSTFILEnc="${NAME}_${im1}.restart.nc"
    fi   
    
    
    echo ''
    echo "start: ${NAME}_elmer_$i from $RSTFILEnc"
    echo '======' 
    echo ''
    
    ## simulation with or without smb anomaly
    SMB_lc=$(echo "$SMB_METHOD" | tr '[:upper:]' '[:lower:]')
    case "$SMB_lc" in
             anomaly)
		    echo "Simulation with constant SMB and variable ASMB" 
                    Ini_Condition="! asmb = Real 0.0"
		    Execution1="Before simulation"
		    Execution2="Before Timestep"
		    DATA="asmb"
                    ;;
             constant)
                    echo "Simulation with constant SMB" 
                    Ini_Condition="asmb = Real 0.0"
                    Execution1="Before simulation"
		    Execution2="Never"
		    DATA="asmb"
                    ;;
             variable)
                    echo "Simulation with variable SMB" 
                    Ini_Condition="asmb = Real 0.0"
		    Execution1="Never"
                    Execution2="Before Timestep"
                    DATA="smb"
                    ;;
             *)     echo "Sorry SMB method not found; exit"; exit 42
    esac
 
    ## parameterisation of basal melting
    PARAM_MELT_lc=$(echo "$MELT" | tr '[:upper:]' '[:lower:]')
    case "$PARAM_MELT_lc" in
             pico)
                    echo "Simulation with PICO parameterisation"
		    # sed -e '/<PARAMETER_MELT>/ {' -e 'r PICO_Param.txt' -e 'd' -e '}' -i $WELMER/elmer_t${i}.sif
		    sed '/<PARAMETER_MELT>/ s/.*/cat PICO_Param.txt/e' ${NAME}_elmer.sif > ${NAME}_elmer_melt.sif
                    ;;
             quadratic)
                    echo "Simulation with QUADRATIC parameterisation"
		    sed '/<PARAMETER_MELT>/ s/.*/cat QUADRATIC_Param.txt/e' ${NAME}_elmer.sif > ${NAME}_elmer_melt1.sif
		    sed -e "s/<SLOPE>/$SLOPE/g" \
		        -e "s/<CORR_T>/$CORRECTION/g" ${NAME}_elmer_melt1.sif > ${NAME}_elmer_melt.sif
		    rm -f ${NAME}_elmer_melt1.sif
		    ;;
             *)     echo "Sorry Melt Parameterisation not found; exit"; exit 42
    esac



    # prepare sif
    sed -e "s/<ID-1>/${im1}/g"               \
        -e "s/<ID>/${i}/g"                   \
        -e "s/<NSTEPS>/$NSTEP/g"             \
        -e "s/<STPINDAYS>/$TIME_STP/g"       \
        -e "s/<STARTYEAR>/$START_SIMU/g"     \
	-e "s/<MELT>/$MELT/g"                \
        -e "s/<OFFSET>/$OFFSET/g"            \
        -e "s/<OFFSETOC>/$OFFSETOC/g"        \
	-e "s/<INI_COND>/${Ini_Condition}/g" \
	-e "s/<EXEC_C>/${Execution1}/g"      \
	-e "s/<EXEC_V>/${Execution2}/g"      \
	-e "s/<DATA_V>/${DATA}/g"            \
        -e "s/<RSTFILEnc>/$RSTFILEnc/g"      ${NAME}_elmer_melt.sif  > $WELMER/elmer_t${i}.sif  

    # prepare run script
    sed -e "s!<NAME>!${NAME}_$i!g"       \
        -e "s!<NNODES>!${NN}!g"          \
        -e "s!<GROUPUSR>!${GROUPUSR}!g"  \
        -e "s!<WALLTIME>!${WALLTIME}!g"  \
        -e "s!<NTASKS>!${NP}!g"          run_arch.slurm > run_elmer_${i}.slurm

    sed -e "s!<RSTFILEnc>!$RSTFILEnc!g"  \
        -e "s!<ECONFIG>!$CONFIG!g"       \
        -e "s!<ECASE>!$CASE!g"           \
        -e "s!<HELMER>!$HELMER!g"        \
        -e "s!<ID>!${i}!g"               \
        -e "s!<ID-1>!${im1}!g"           \
        -e "s!<ID+1>!${ip1}!g"           run_elmer_skel.bash >> run_elmer_${i}.slurm
   




    ## create Friction_Material
    FRICTION_lc=$(echo "$FRICTION" | tr '[:upper:]' '[:lower:]')
    case "$FRICTION_lc" in
             linear)
                    echo "SSA Friction Law = String \"linear\"" > $WELMER/Friction_Material.sif || exit 42
                    echo "SSA Friction Parameter = Equals frC" >> $WELMER/Friction_Material.sif || exit 42
                    XIOS_frc_def="\<field id=\"frc\"  long_name=\"ssa_linear_friction_coefficent\" unit=\"MPa m-1 d\"    grid_ref=\"GridNodes\" \/\>"
                    ;;
             weertman)
                    echo "SSA Friction Law = String \"Weertman\"" > $WELMER/Friction_Material.sif || exit 42
                    echo "SSA Friction Parameter = Equals frC" >> $WELMER/Friction_Material.sif || exit 42
                    echo "SSA Friction Exponent = Real $ 1.0/n" >> $WELMER/Friction_Material.sif || exit 42
                    echo "SSA Friction Linear Velocity = Real $ Vmin" >> $WELMER/Friction_Material.sif || exit 42
                    n=$(\echo $(\grep "^ *\$n *= *[0-9]" ANT50.GL1-ISMIP6_elmer.param | \cut -d= -f2 | \cut -d! -f1) | \sed s/.0$//g)
		    echo "n = $n"
		    XIOS_frc_def="\<field id=\"frc\"  long_name=\"ssa_weertman_friction_coefficent\" unit=\"MPa m^-1/$n d^1/$n\"    grid_ref=\"GridNodes\" \/\>"
                    ;;
             "regularized coulomb")
                    echo "SSA Friction Law = String \"Regularized coulomb\"" > $WELMER/Friction_Material.sif || exit 42
                    echo "SSA Friction Parameter = Variable frC, haf" >> $WELMER/Friction_Material.sif || exit 42
                    echo "   Real Procedure "SlipCoef" "Calcul_Slc"" >> $WELMER/Friction_Material.sif || exit 42
                    echo "SSA Friction Exponent = Real $ 1.0/n" >> $WELMER/Friction_Material.sif || exit 42
                    echo "SSA Friction Linear Velocity = Real $ Vmin" >> $WELMER/Friction_Material.sif || exit 42
                    echo "SSA Friction Threshold Velocity=Real $ Vthres" >> $WELMER/Friction_Material.sif  || exit 42
                    echo "SSA Friction Threshold Height=Real $ Hthres" >> $WELMER/Friction_Material.sif  || exit 42
                    XIOS_frc_def="\<field id=\"frc\"  long_name=\"ssa_RegCoulomb_friction_coefficent\" unit=\"MPa\"    grid_ref=\"GridNodes\" \/\>"
                    ;;
 
              *)  echo "Sorry Friction Law not ready; exit"; exit 42
    esac

    echo "----------------------------------------------------"
    echo "Friction law definition :"
    cat $WELMER/Friction_Material.sif
    echo "----------------------------------------------------"


    # prepare xml file
    year=`grep '$yearinsec' ${CONFIG}-${CASE}_elmer.param`
    year=${year#*=}      #supress before =
    year=${year%\!*}     #supress after !
    
    day=`grep '$dayinsec' ${CONFIG}-${CASE}_elmer.param`
    day=${day#*=}      #supress before =
    day=${day%\!*}     #supress after !

    rho=`grep '$rhoi_SI' ${CONFIG}-${CASE}_elmer.param`
    rho=${rho#*=}      #supress before =
    rho=${rho%\!*}     #supress after !

    sed -e "s!<SEC_YEAR>!${year}!g"                \
        -e "s!<SEC_DAY>!${day}!g"                  \
     	-e "s!<RHOI>!$rho!g"                       \
        -e "s!<ORIGINE>!$TIME_ORIGINE!g"           \
	-e "s!<FRICTION_FIELD>!${XIOS_frc_def}!g"  context_elmer.xml > $WELMER/context_elmer.xml || exit 42
    #sed "s@<FRICTION_FIELD>@$XIOS_frc_def@g" context_elmer.xml > $WELMER/context_elmer.xml || exit 42
    
    # manage status file
    if [[ $((ijob)) == 1 ]];then touch ${HELMER}/zELMER_${i}_READY; fi

    # submit job
    if [ ! -z "$jobid0" ];then
        jobid=$(submit_elmer_dependency $jobid0)
        echo "        id        : $jobid"
        echo "        dependency: $jobid0"
    else
        jobid=$(submit_elmer)
        echo "        id        : $jobid"
    fi
    jobid0=$jobid
 
done
rm -f ${NAME}_elmer_melt.sif
echo ''
