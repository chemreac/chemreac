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
  CXXFLAGS=-std=c++14 -DDEBUG
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
  CXXFLAGS=-std=c++14
  LIBS=-lrt
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
FLAGS += -D_GLIBCXX_DEBUG

TARGETS=test_chemreac.out test_chemreac_sundials.out test_chemreac_catch.out
GENERATED=chemreac.cpp test_chemreac test_chemreac.out test_chemreac_sundials test_chemreac_sundials.out test_utils.o chemreac.o

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

PYTHON?=python
FLAGS+=$(EXTRA_COMPILE_ARGS)
LIBS+=$(shell $(PYTHON) -c "from pycvodes import config as c; print(' '.join(['-l%s' % l for l in (c['LAPACK'].split(',') + c['SUNDIALS_LIBS'].split(','))]))")

LOCALLIB=/usr/local/lib
LIBDIRS=-L${LOCALLIB}
#note if missing libraries do, e.g.:
#env LD_LIBRARY_PATH=/usr/lib/gcc/x86_64-linux-gnu/4.6/:$LD_LIBRARY_PATH make

INCLUDE ?= -I../chemreac/include -I$(shell $(PYTHON) -c "import finitediff as fd; print(fd.get_include())") -I$(shell $(PYTHON) -c "import block_diag_ilu as bdi; print(bdi.get_include())") -I$(shell $(PYTHON) -c "import pycvodes as pc; print(pc.get_include())")

all: $(TARGETS)

.PHONY: clean all

clean:
	rm $(GENERATED)

chemreac.cpp: ../src/chemreac_template.cpp
	$(PYTHON) enmako.py -g "{\"WITH_OPENMP\": $(WITH_OPENMP_TRFL)}" -o $@ $<


test_chemreac: chemreac.cpp test_chemreac.cpp test_utils.o
	$(CXX) $(CXXFLAGS) $(FLAGS) $(INCLUDE) -o $@ $^ $(LIBS) # -lrt only needed if we're not using OpenMP


test_chemreac_sundials: chemreac.o test_utils.o test_chemreac_sundials.cpp
	$(CXX) $(CXXFLAGS) $(FLAGS) $(INCLUDE) $(LIBDIRS) -o $@ $^ -lm $(LIBS)

gmon.out: test_chemreac
	./$<

profile_test_chemreac.out: test_chemreac gmon.out
	./$<
	gprof ./$< | tee $@

test_chemreac.out: test_chemreac
	$(CONTEXT) ./$< | tee $@

test_chemreac_sundials.out: test_chemreac_sundials
	$(CONTEXT) ./$< >$@

test_chemreac_catch.out: test_chemreac_catch
	$(CONTEXT) ./$< | tee $@

catch.hpp: catch.hpp.bz2
	bunzip2 -k -f $<

test_chemreac_catch: test_chemreac_catch.cpp catch.hpp chemreac.o
	$(CXX) $(CXXFLAGS) $(FLAGS) $(INCLUDE) $(LIBDIRS) -o $@ $< chemreac.o test_utils.o $(LIBS)

%.o: %.cpp
	$(CXX) -fPIC $(CXXFLAGS) $(FLAGS) $(INCLUDE) -c -o $@ $<

