include config.mk
FC = g++

FFLAGS = -O3 -I$(OPT) -I/usr/include/suitesparse
OTHER = -L$(LIB) -lmfem -lnlopt -lm

all: surface
surface: minimalsurface-example.cpp minimalsurface.cpp
	$(FC) $(FFLAGS) $^ -o $@ $(OTHER)	

clean:
	rm -f surface *.o



