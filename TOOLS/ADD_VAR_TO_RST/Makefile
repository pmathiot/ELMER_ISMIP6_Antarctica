
SRC=./

BLD=./

NCINC=-I$(NETCDFFORTRAN_INCDIR) -I$(NETCDFC_INCDIR)

NCLIB=-L$(NETCDFFORTRAN_LIBDIR) -L$(NETCDFC_LIBDIR) -lnetcdf -lnetcdff

VPATH=$(BLD)

FFLAGS=-I$(XIOS_INC) -g -traceback -O0 -fpe0 -CB -ftrapuv -debug all -CA
# -fpe0 -CB -ftrapuv -debug

LIBS= USF_SSAViscosity nearestpoint

all: $(LIBS)

#---------------------------------------------------
# USR function
USF_SSAViscosity:$(SRC)/USF_SSAViscosity.F90
	elmerf90 $^ -o $(BLD)/$@ $(FFLAGS) $(ELMERINC)


#---------------------------------------------------
# PP
# PP define basin
nearestpoint:$(SRC)/nearestpoint.F90
	elmerf90 $(NCINC) $^ -o $(BLD)/$@ $(FFLAGS) $(NCLIB)

clean:
	cd $(BLD) && rm -f $(LIBS)
