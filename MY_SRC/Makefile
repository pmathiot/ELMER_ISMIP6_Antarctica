
SRC=./

BLD=../MY_BLD/

NCINC=-I$(NETCDFFORTRAN_INCDIR) -I$(NETCDFC_INCDIR)

NCLIB=-L$(NETCDFFORTRAN_LIBDIR) -L$(NETCDFC_LIBDIR) -lnetcdf -lnetcdff

VPATH=$(BLD)

FFLAGS=-I$(XIOS_INC) -g -traceback -O0 -fpe0 -CB -ftrapuv -debug all -CA
# -fpe0 -CB -ftrapuv -debug

LIBS= USF_WaterPressure USF_SSAViscosity USF_ExitCondition NCFILE_BASALMELT nearestpoint INITMIP_Scalar_OUTPUT boxmodel_solver FrontThickness_mask USF_CondFront USF_CondDistance SlipCoef
#LIBS= USFs_Ronnie_inversion INITMIP_Scalar_OUTPUT Save_Netcdf Flotation3 Flotation2 boxmodel_solver nearestpoint FrontThickness_mask USF_CondFront USF_CondDistance USF_WaterPressure PDC_BASALMELT Elmer2NEMO_Interp NCFILE_BASALMELT TransientRestartInit USF_ExitCondition SlipCoef

all: $(LIBS)

#---------------------------------------------------
# USR function
USF_WaterPressure:$(SRC)/USF_WaterPressure.F90
	elmerf90 $^ -o $(BLD)/$@ $(FFLAGS)

USF_ExitCondition:${SRC}/USF_ExitCondition.F90
	elmerf90 $^ -o $(BLD)/$@ $(FFLAGS)

USF_SSAViscosity:$(SRC)/USF_SSAViscosity.F90
	elmerf90 $^ -o $(BLD)/$@ $(FFLAGS) $(ELMERINC)


#---------------------------------------------------
# PP
# PP define basin
nearestpoint:$(SRC)/nearestpoint.F90
	elmerf90 $(NCINC) $^ -o $(BLD)/$@ $(FFLAGS) $(NCLIB)

# PP compute scalar per basin
INITMIP_Scalar_OUTPUT:$(SRC)/INITMIP_Scalar_OUTPUT_modif.F90
	elmerf90 $(NCINC) $^ -o $(BLD)/$@ $(FFLAGS) $(NCLIB)


#---------------------------------------------------
# BMB 
# BMB: PICO
boxmodel_solver:$(SRC)/boxmodel_solver.F90
	elmerf90 $(NCINC) $^ -o $(BLD)/$@ $(FFLAGS) $(NCLIB)

FrontThickness_mask:$(SRC)/FrontThickness_mask.F90
	elmerf90 $^ -o $(BLD)/$@ $(FFLAGS)

USF_CondFront:$(SRC)/USF_CondFront.F90
	elmerf90 $^ -o $(BLD)/$@ $(FFLAGS)

USF_CondDistance:$(SRC)/USF_CondDistance.F90
	elmerf90 $^ -o $(BLD)/$@ $(FFLAGS)

# BMB: read a NEMO melt in netcdf
NCFILE_BASALMELT: $(SRC)/NCFILE_BASALMELT.F90
	elmerf90 $(NCINC) $^ -o $(BLD)/$@ $(FFLAGS) $(NCLIB)

#---------------------------------------------------
# Slip Coeficient correction Coulomb
SlipCoef:$(SRC)/SlipCoef.F90
	elmerf90 $^ -o $(BLD)/$@ $(FFLAGS)

clean:
	cd $(BLD) && rm -f $(LIBS)
