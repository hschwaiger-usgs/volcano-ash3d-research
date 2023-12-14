###############################################################################
##########  Intel Fortran Compiler  ###########################################
SYSFPPFLAG = -DIFORT

FCHOME = /opt/intel/oneapi/compiler/latest/linux/
FC = $(FCHOME)/bin/intel64/ifort

#LAPACKLIB = -L$(FCHOME)/lib64 -llapack
#LAPACKLIB = -L$(HOME)/Programs/Tarballs/lapack-3.11 -llapack
LAPACKLIB = -L/opt/intel/oneapi/mkl/latest/lib/intel64 -qmkl

COMPINC = -I./ -I$(FCHOME)/include
COMPLIBS = -L./ -L$(FCHOME)/lib -L/opt/intel/oneapi/mkl/latest/lib/intel64 -lirc -limf
LIBS = $(LAPACKLIB) $(COMPLIBS) $(USGSLIBDIR) $(USGSINC) $(COMPINC) $(USGSLIB) $(DATALIBS)

# Debugging flags
ifeq ($(RUN), DEBUG)
    FFLAGS = -debug full -g3 -O0 -CB -init:snan,arrays -warn all -gen-interfaces -traceback -check all -check bounds -ftrapuv
endif
ifeq ($(RUN), DEBUGOMP)
    FFLAGS = -debug full -g3 -O0 -CB -init:snan,arrays -warn all -gen-interfaces -traceback -check all -check bounds -ftrapuv -qopenmp -parallel
endif
# Profiling flags
ifeq ($(RUN), PROF)
    FFLAGS = -g2 -pg
endif
# Production run flags
ifeq ($(RUN), OPT)
    FFLAGS = -O3 -unroll -unroll-aggressive -ftz -inline-forceinline -finline-functions -ipo -ip
endif
ifeq ($(RUN), OMPOPT)
    FFLAGS = -O3 -unroll -unroll-aggressive -ftz -inline-forceinline -finline-functions -ipo -ip -qopenmp -parallel
endif
  # Preprocessing flags
FPPFLAGS =  -fpp $(OSFPPFLAG) $(DATAFPPFLAG) $(LIMFPPFLAG) $(DIFFFPPFLAG) $(memFPPFLAG) $(extdatFPPFLAG) $(FASTFPPFLAG) $(PIIFPPFLAG) $(SYSFPPFLAG) $(PLOTFPPFLAG)
  # Extra flags
EXFLAGS =
###############################################################################
