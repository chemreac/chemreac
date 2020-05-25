#
# Usage example:
#
#  $ make -B CONTEXT="gdb -ex run --args" EXTRA_COMPILE_ARGS='-D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC'
#

CXX?=g++ 
CC?=gcc
PERFFLAGS=-funroll-loops -O3 -ffast-math # -march=native
PROFFLAGS=-pg
CXXFLAGS=-std=c++11 -Wall -Wextra -pedantic -Werror $(EXTRA_COMPILE_ARGS) ${CFLAGS}
OPENMPLIBS=-lgomp
OPENMPFLAG=-fopenmp
LIBS=-lrt -llapack -lblas
SUNDIALS_LIBS ?= $(${PYTHON} -c "from pycvodes._libs import print_libs_linkline as pll; pll()")
ifeq ($(OPTIMIZE),1)
  CONTEXT ?= # /usr/bin/time
  CXXFLAGS += $(PERFFLAGS) -DBENCHMARK
else
  CXXFLAGS +=-g -DDEBUG 
endif

GENERATED=chemreac.cpp test_chemreac test_utils.o chemreac.o

ifeq ($(PROFILING),1)
  CXXFLAGS += $(PROFFLAGS) -DPROFILING 
  GENERATED += profile_test_chemreac.out gmon.out
endif
ifeq ($(WITH_OPENMP),1)
  CXXFLAGS += $(OPENMPFLAG)
  LIBS += $(OPENMPLIBS)
  WITH_OPENMP_TRFL="True"
else
  WITH_OPENMP_TRFL="False"
endif

CXXFLAGS += $(EXTRA_FLAGS)


INCLUDE ?= -I../chemreac/include -I../external/anyode/include -I$(shell python3 -c "import finitediff as fd; print(fd.get_include())") -I$(shell python3 -c "import block_diag_ilu as bdi; print(bdi.get_include())") -I$(shell python3 -c "import pycvodes as pc; print(pc.get_include())")

%.o: %.cpp
	$(CXX) -fPIC $(CXXFLAGS) $(INCLUDE) -c -o $@ $<


.PHONY: all test clean

all: test $(GENERATED)

test: test_chemreac test_chemreac_sundials test_chemreac_catch
	$(CONTEXT) ./test_chemreac
	$(CONTEXT) ./test_chemreac_sundials
	./test_chemreac_catch

clean:
	rm $(GENERATED)

chemreac.cpp: ../src/chemreac.cpp.mako
	python3 enmako.py -g "{\"WITH_OPENMP\": $(WITH_OPENMP_TRFL)}" -o $@ $<


test_chemreac: chemreac.cpp test_chemreac.cpp test_utils.o
	$(CXX) $(CXXFLAGS) $(INCLUDE) \
	 -o $@ $^ $(LIBS) # -lrt only needed if we're not using OpenMP


test_chemreac_sundials: chemreac.o test_utils.o test_chemreac_sundials.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $^ $(LDFLAGS) $(SUNDIALS_LIBS) -lm

gmon.out: test_chemreac
	./$<

profile_test_chemreac.out: test_chemreac gmon.out
	./$<
	gprof ./$< > $@

catch.hpp: catch.hpp.bz2
	bzcat $< >$@

test_chemreac_catch: test_chemreac_catch.cpp catch.hpp test_utils.o chemreac.o
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $< chemreac.o test_utils.o $(LIBS)
