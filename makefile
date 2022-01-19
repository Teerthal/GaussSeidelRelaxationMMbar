FC=gfortran
#FSTRICT=-Wextra -Wno-unused -Wsurprising -Wall #-Wimplicit #-pedantic
#FVECTORIZE=-ftree-vectorizer-verbose=2 -ftree-vectorize -ftree-vect-loop-version
#FOPTIMIZE=-fvariable-expansion-in-unroller -funroll-loops $(FVECTORIZE)
#OPTS = -O3 -g $(FOPTIMIZE) #-fopenmp


FFLAGS = $(OPTS) $(PROF)
LDFLAGS = $(FFLAGS)

# 
OBJ = relaxMain.o \
      derivatives6thOrder.o \
      relaxEuler.o \
      relaxEnergySO3.o \
      covariantDerivsSO3.o \
      fieldStrengthsSO3.o \
      fluxesNumRelSO3.o \
      coeff.o \
      icMMbarTwistedNumerical.o \
#      icMonopole.o \
#      velfactors.o \
#      energySO3.o \
#      relaxChangeName.o \
#      helicitySO3.o \
# crankNicholsonLeapforward.o \
#      crankNicholsonMain.o \
#      evolveeuler.o \
# relaxLeapForward.o \

all: cnSO3.out 

.f.o:
	$(FC) $(FFLAGS) -c $*.f

cnSO3.out: $(OBJ) 
	$(FC) -o  $@ $(LDFLAGS) $(OBJ) 

cleanf:
	rm: -f *.o *.dvi *.aux *.log core.* *~
