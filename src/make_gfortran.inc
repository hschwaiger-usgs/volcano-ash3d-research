###############################################################################
##########  GNU Fortran Compiler  #############################################
SYSFPPFLAG = -DGFORTRAN

FCHOME = /usr
FC = $(FCHOME)/bin/gfortran

#LAPACKLIB = -L$(FCHOME)/lib64 -llapack
LAPACKLIB = -L$(HOME)/Programs/Tarballs/lapack-3.11 -llapack

COMPINC = -I./ -I$(FCHOME)/lib64/gfortran/modules
COMPLIBS = -L$(FCHOME)/lib64
LIBS = $(LAPACKLIB) $(COMPLIBS) $(USGSLIBDIR) $(USGSINC) $(COMPINC) $(USGSLIB) $(DATALIBS)

# Debugging flags
ifeq ($(RUN), DEBUG)
    # consider including -lefence for debugging
    FFLAGS = -O0 -g3 -Wall -Wextra -fimplicit-none -Wline-truncation -Wcharacter-truncation -Waliasing -Wimplicit-interface -Wunused-parameter -fwhole-file -fcheck=all -pedantic -fbacktrace -Wunderflow -ffpe-trap=invalid,zero,overflow -Wuninitialized -Wmaybe-uninitialized -finit-integer=-2147483647 -finit-real=snan -fno-unsafe-math-optimizations -fsignaling-nans -frounding-math -Wno-surprising
endif
ifeq ($(RUN), DEBUGOMP)
    FFLAGS = -O0 -g3 -Wall -Wextra -fimplicit-none -Wline-truncation -Wcharacter-truncation -Waliasing -Wimplicit-interface -Wunused-parameter -fwhole-file -fcheck=all -pedantic -fbacktrace -Wunderflow -ffpe-trap=invalid,zero,overflow -Wuninitialized -Wmaybe-uninitialized -finit-integer=-2147483647 -finit-real=snan -fno-unsafe-math-optimizations -fsignaling-nans -frounding-math -Wno-surprising -ffpe-summary=all -std=f2008 -fopenmp -lgomp
endif
# Profiling flags
ifeq ($(RUN), PROF)
    FFLAGS = -g -pg -w -fno-math-errno -funsafe-math-optimizations -fno-trapping-math -fno-signaling-nans -fcx-limited-range -fno-rounding-math
endif
# Production run flags
ifeq ($(RUN), OPT)
    FFLAGS = -O3 -w -fno-math-errno -funsafe-math-optimizations -fno-trapping-math -fno-signaling-nans -fcx-limited-range -fno-rounding-math -ftree-vectorize -funroll-loops -finline-functions
endif
ifeq ($(RUN), OMPOPT)
    FFLAGS = -O3 -w -fno-math-errno -funsafe-math-optimizations -fno-trapping-math -fno-signaling-nans -fcx-limited-range -fno-rounding-math -ftree-vectorize -funroll-loops -finline-functions -fopenmp -lgomp
endif
  # Preprocessing flags
FPPFLAGS = -x f95-cpp-input $(OSFPPFLAG) $(DATAFPPFLAG) $(LIMFPPFLAG) $(DIFFFPPFLAG) $(memFPPFLAG) $(extdatFPPFLAG) $(FASTFPPFLAG) $(PIIFPPFLAG) $(SYSFPPFLAG) $(PLOTFPPFLAG)
  # Extra flags
EXFLAGS =
###############################################################################
