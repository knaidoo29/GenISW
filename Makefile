CC = g++
OPTIMIZE = -c -O3

GSL_INC = -I/Users/krishna/Programs/build/gsl/include
GSL_LIB = #-L/Users/krishna/Programs/build/gsl/lib -lgsl
BOOST_INC = -I/Users/krishna/Programs/build/boost/include
BOOST_LIB = #-L/Users/krishna/Programs/build/boost/lib

INC = $(GSL_INC) $(BOOST_INC)
LIBS = $(GSL_LIB) $(BOOST_LIB)
CFLAGS = $(OPTIMIZE)

UTIL = src/util_interp.o src/util_omp.o src/util_param.o src/util_progress.o src/util_read.o
OBJS = src/cosmo.o src/spherical_bessel.o $(UTIL)

all: lib/libgenisw.a

lib/libgenisw.a: $(OBJS)
	 ar cr lib/libgenisw.a $(OBJS)

.cpp.o:
	 $(CC) $(CFLAGS) $(INC) $< $(LIBS) -o $@

clean:
	 rm src/*.o*
	 rm lib/*.a*
	 rm .DS_Store
