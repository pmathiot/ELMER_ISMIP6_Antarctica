
SRC=./

BLD=../MY_BLD/

NCINC=-I$(NETCDFFORTRAN_INCDIR) -I$(NETCDFC_INCDIR)

NCLIB=-L$(NETCDFFORTRAN_LIBDIR) -L$(NETCDFC_LIBDIR) -lnetcdf -lnetcdff

VPATH=$(BLD)

FFLAGS=-I$(XIOS_INC) -O2 -fpe0 -ftrapuv -g -traceback

LIBS= USF_ExitCondition boxmodel_solver_t FrontThickness_mask USF_CondFront USF_CondDistance SlipCoef SSASolver_bugfix UGridDataReader_t nearestpoint

all: $(LIBS)

#---------------------------------------------------
# USR function
USF_ExitCondition:${SRC}/USF_ExitCondition.F90
	elmerf90 $^ -o $(BLD)/$@ $(FFLAGS)

#---------------------------------------------------
# Post Porcessing
# PP define basin
nearestpoint:$(SRC)/nearestpoint.F90
	elmerf90 $(NCINC) $^ -o $(BLD)/$@ $(FFLAGS) $(NCLIB)
#
#---------------------------------------------------

# BMB 
# BMB: PICO
boxmodel_solver_t:$(SRC)/boxmodel_solver_t.F90
	elmerf90 $(NCINC) $^ -o $(BLD)/$@ $(FFLAGS) $(NCLIB)

FrontThickness_mask:$(SRC)/FrontThickness_mask.F90
	elmerf90 $^ -o $(BLD)/$@ $(FFLAGS)

USF_CondFront:$(SRC)/USF_CondFront.F90
	elmerf90 $^ -o $(BLD)/$@ $(FFLAGS)

USF_CondDistance:$(SRC)/USF_CondDistance.F90
	elmerf90 $^ -o $(BLD)/$@ $(FFLAGS)

#---------------------------------------------------
#SMB
UGridDataReader_t:$(SRC)/UGridDataReader_t.F90
	elmerf90 $(NCINC) $^ -o $(BLD)/$@ $(FFLAGS) $(NCLIB)

#---------------------------------------------------
# Slip Coeficient correction Coulomb
SlipCoef:$(SRC)/SlipCoef.F90
	elmerf90 $^ -o $(BLD)/$@ $(FFLAGS)

clean:
	cd $(BLD) && rm -f $(LIBS)

#---------------------------------------------------
# Bug fix SSA coulomb reg
SSASolver_bugfix:$(SRC)/SSASolver_bugfix.F90
	elmerf90 $^ -o $(BLD)/$@ $(FFLAGS)
clean:
	cd $(BLD) && rm -f $(LIBS)
