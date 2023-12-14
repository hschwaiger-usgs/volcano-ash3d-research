###############################################################################
##########  AMD Optimizing C/C++/Fortran Compiler (aocc) ######################
SYSFPPFLAG = -DAOCC

FCHOME = /opt/AMD/aocc-compiler-4.0.0
FC = $(FCHOME)/bin/flang

LAPACKLIB = -L$(FCHOME)/lib64 -llapack
#LAPACKLIB = -L$(HOME)/Programs/Tarballs/lapack-3.11 -llapack

COMPINC = -I./ -I$(FCHOME)/include
COMPLIBS = -L./ -L$(FCHOME)/lib
LIBS = $(LAPACKLIB) $(COMPLIBS) $(USGSLIBDIR) $(USGSINC) $(COMPINC) $(USGSLIB) $(DATALIBS)

# Debugging flags
ifeq ($(RUN), DEBUG)
    FFLAGS = -O0 -g3 -Wall -Mstandard -fPIC -v
endif
ifeq ($(RUN), DEBUGOMP)
    FFLAGS = -O0 -g3 -Wall -Mstandard -fPIC -v -fopenmp
endif
# Profiling flags
ifeq ($(RUN), PROF)
    FFLAGS = -g2 -pg -Mstandard -fPIC
endif
# Production run flags
ifeq ($(RUN), OPT)
    #FFLAGS = -O3 -zopt -ftz -w -ipo
    FFLAGS = -Ofast -Mstandard -fPIC -mavx -fplugin-arg-dragonegg-llvm-codegen-optimize=3 -fplugin-arg-dragonegg-llvm-ir-optimize=3
endif
ifeq ($(RUN), OMPOPT)
    FFLAGS = -Ofast -Mstandard -fPIC -mavx -fplugin-arg-dragonegg-llvm-codegen-optimize=3 -fplugin-arg-dragonegg-llvm-ir-optimize=3 -fopenmp
endif
  # Preprocessing flags
FPPFLAGS =  -Mpreprocess $(OSFPPFLAG) $(DATAFPPFLAG) $(LIMFPPFLAG) $(DIFFFPPFLAG) $(memFPPFLAG) $(extdatFPPFLAG) $(FASTFPPFLAG) $(PIIFPPFLAG) $(SYSFPPFLAG) $(PLOTFPPFLAG)
  # Extra flags
EXFLAGS =
###############################################################################
