ifeq ($(COMPILER_VENDOR),intel) # 2013 SP 1
  FC=ifort
  CXX=icpc 
  CC=icc
  FORTCPPLINKER=ifort
  FORTCPPLIBS=
  FORTCPPFLAGS=-nofor_main -cxxlib
  FORTFLAGS=-stand f08 -warn all
  FLAGS=-Wall
  PROFFLAGS=-prof-gen -prof-dir=./profdata
  PERFFLAGS=-fast
  CXXFLAGS=-std=c++11 -DDEBUG
  OPENMPLIBS=-liomp5
  OPENMPFLAG=-openmp
  C99=-std=c99
else
  FC=gfortran
  CXX=g++ 
  CC=gcc
  FORTCPPLINKER=g++
  FORTCPPLIBS=-lgfortranbegin -lgfortran
  FORTCPPFLAGS=
  FORTFLAGS=-std=f2008 
  PERFFLAGS=-funroll-loops -O3 -ffast-math # -march=native
  PROFFLAGS=-pg
  FLAGS=-Wall -Wextra -pedantic
  OPENMPLIBS=-lgomp
  OPENMPFLAG=-fopenmp
  CXXFLAGS=-std=c++11
  LIBS=-lrt -llapack
  C99=-std=c99
endif

ifeq ($(OPTIMIZE),1)
  CONTEXT = time
  FLAGS += $(PERFFLAGS) -DBENCHMARK
  FORTFLAGS += $(PERFFLAGS)
else
  CONTEXT = valgrind
  FLAGS +=-g -DDEBUG 
  FORTFLAGS +=-g #-DDEBUG
  ifeq ($(COMPILER_VENDOR),intel) # 2013 SP 1
    FLAGS +=-debug all
    FORTFLAGS +=-debug all
  endif
endif

TARGETS=test_chemreac.out test_chemreac_sundials.out test_chemreac_catch.out
GENERATED=chemreac.cpp test_chemreac test_chemreac.out fornberg.o test_chemreac_sundials test_chemreac_sundials.out test_utils.o chemreac.o c_fornberg.o fornberg.mod c_fornberg.mod

ifeq ($(PROFILING),1)
  FLAGS += $(PROFFLAGS) -DPROFILING 
  FORTFLAGS += $(PROFFLAGS) -DPROFILING 
  TARGETS += profile_test_chemreac.out
  GENERATED += profile_test_chemreac.out gmon.out
endif
ifeq ($(WITH_OPENMP),1)
  FLAGS += $(OPENMPFLAG)
  LIBS += $(OPENMPLIBS)
  WITH_OPENMP_TRFL="True"
else
  WITH_OPENMP_TRFL="False"
endif

FLAGS += $(EXTRA_FLAGS)


LOCALLIB=/usr/local/lib
LIBDIRS=-L${LOCALLIB}

#note if missing libraries do, e.g.:
#env LD_LIBRARY_PATH=/usr/lib/gcc/x86_64-linux-gnu/4.6/:$LD_LIBRARY_PATH make

INCLUDE ?= -I.. -I../finitediff/finitediff

%.o: %.cpp
	$(CXX) -fPIC $(CXXFLAGS) $(FLAGS) $(INCLUDE) -c -o $@ $<

all: $(TARGETS)

.PHONY: clean all

clean:
	rm $(GENERATED)

chemreac.cpp: ../chemreac_template.cpp
	python enmako.py -g "{\"WITH_OPENMP\": $(WITH_OPENMP_TRFL)}" -o $@ $<

fornberg.o: ../finitediff/finitediff/fornberg.f90
	$(FC) $(FORTFLAGS) -c -o $@ $^

c_fornberg.o: ../finitediff/finitediff/c_fornberg.f90
	$(FC) $(FORTFLAGS) -c -o $@ $^

test_chemreac: chemreac.cpp test_chemreac.cpp fornberg.o c_fornberg.o test_utils.o
	$(CXX) $(CXXFLAGS) $(FLAGS) $(INCLUDE) \
	 -o $@ $^ $(LIBS) $(FORTCPPLIBS) # -lrt only needed if we're not using OpenMP


test_chemreac_sundials: chemreac.o fornberg.o c_fornberg.o test_utils.o test_chemreac_sundials.cpp
	$(CXX) $(CXXFLAGS) $(FLAGS) $(INCLUDE) $(LIBDIRS) -o $@ $^ -lm $(FORTCPPLIBS) -lsundials_cvodes -llapack -lsundials_nvecserial

gmon.out: test_chemreac
	./$<

profile_test_chemreac.out: test_chemreac gmon.out
	./$<
	gprof ./$< > $@

test_chemreac.out: test_chemreac
	$(CONTEXT) ./$< >$@

test_chemreac_sundials.out: test_chemreac_sundials
	$(CONTEXT) ./$< >$@

test_chemreac_catch.out: test_chemreac_sundials
	$(CONTEXT) ./$< >$@

catch.hpp: catch.hpp.bz2
	bunzip2 -k -f $<

test_chemreac_catch: test_chemreac_catch.cpp catch.hpp
	$(CXX) $(CXXFLAGS) $(FLAGS) $(INCLUDE) $(LIBDIRS) -o $@ $< chemreac.o fornberg.o c_fornberg.o test_utils.o $(LIBS) $(FORTCPPLIBS)
