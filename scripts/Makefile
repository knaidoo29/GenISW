CC = g++
CCOMP = /Users/krishna/Programs/build/openmpi/bin/mpic++

OPTIMIZE = -O3

GENISW_INC = -I/Users/krishna/Research/github/GenISW/include
GENISW_LIB = -L/Users/krishna/Research/github/GenISW/lib -lgenisw
GSL_INC = -I/Users/krishna/Programs/build/gsl/include
GSL_LIB = -L/Users/krishna/Programs/build/gsl/lib -lgsl
BOOST_INC = -I/Users/krishna/Programs/build/boost/include
BOOST_LIB = -L/Users/krishna/Programs/build/boost/lib

INC = $(GENISW_INC) $(GSL_INC) $(BOOST_INC)
LIBS = -lm $(GENISW_LIB) $(GSL_LIB) $(BOOST_LIB)
CFLAGS = $(OPTIMIZE)

all: SH2SBT SH2SBT_MPI SBT2ISW SBT2ISW_MPI

SH2SBT:
	$(CC) $(CFLAGS) $(INC) SH2SBT.cpp $(LIBS) -o SH2SBT
	mv SH2SBT ../SH2SBT

SH2SBT_MPI:
	$(CCOMP) $(CFLAGS) $(INC) SH2SBT_MPI.cpp $(LIBS) -o SH2SBT_MPI
	mv SH2SBT_MPI ../SH2SBT_MPI

SBT2ISW:
	$(CC) $(CFLAGS) $(INC) SBT2ISW.cpp $(LIBS) -o SBT2ISW
	mv SBT2ISW ../SBT2ISW

SBT2ISW_MPI:
	$(CCOMP) $(CFLAGS) $(INC) SBT2ISW_MPI.cpp $(LIBS) -o SBT2ISW_MPI
	mv SBT2ISW_MPI ../SBT2ISW_MPI

clean:
	rm -r ../SH2SBT
	rm -r ../SH2SBT_MPI
	rm -r ../SBT2ISW
	rm -r ../SBT2ISW_MPI
