CXX?=g++ 
CC?=gcc
PERFFLAGS=-funroll-loops -O3 -ffast-math # -march=native
PROFFLAGS=-pg
FLAGS=-Wall -Wextra -pedantic -Werror $(EXTRA_COMPILE_ARGS)
OPENMPLIBS=-lgomp
OPENMPFLAG=-fopenmp
CXXFLAGS=-std=c++11
LIBS=-lrt -llapack -lblas

ifeq ($(OPTIMIZE),1)
  CONTEXT ?= time
  FLAGS += $(PERFFLAGS) -DBENCHMARK
else
  FLAGS +=-g -DDEBUG 
endif

TARGETS=test_chemreac.out test_chemreac_sundials.out test_chemreac_catch.out
GENERATED=chemreac.cpp test_chemreac test_chemreac.out test_chemreac_sundials test_chemreac_sundials.out test_utils.o chemreac.o

ifeq ($(PROFILING),1)
  FLAGS += $(PROFFLAGS) -DPROFILING 
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


INCLUDE ?= -I../chemreac/include -I../external/anyode/include -I$(shell python3 -c "import finitediff as fd; print(fd.get_include())") -I$(shell python3 -c "import block_diag_ilu as bdi; print(bdi.get_include())") -I$(shell python3 -c "import pycvodes as pc; print(pc.get_include())")

%.o: %.cpp
	$(CXX) -fPIC $(CXXFLAGS) $(FLAGS) $(INCLUDE) -c -o $@ $<

all: $(TARGETS)

.PHONY: clean all

clean:
	rm $(GENERATED)

chemreac.cpp: ../src/chemreac.cpp.mako
	python3 enmako.py -g "{\"WITH_OPENMP\": $(WITH_OPENMP_TRFL)}" -o $@ $<


test_chemreac: chemreac.cpp test_chemreac.cpp test_utils.o
	$(CXX) $(CXXFLAGS) $(FLAGS) $(INCLUDE) \
	 -o $@ $^ $(LIBS) # -lrt only needed if we're not using OpenMP


test_chemreac_sundials: chemreac.o test_utils.o test_chemreac_sundials.cpp
	$(CXX) $(CXXFLAGS) $(FLAGS) $(INCLUDE) -o $@ $^ $(shell python3 -c 'import pycvodes as p; print(" -l".join([" "] + p.config["SUNDIALS_LIBS"].split(",") + p.config["LAPACK"].split(",")), end="")') -lm

gmon.out: test_chemreac
	./$<

profile_test_chemreac.out: test_chemreac gmon.out
	./$<
	gprof ./$< > $@

test_chemreac.out: test_chemreac
	$(CONTEXT) ./$< >$@

test_chemreac_sundials.out: test_chemreac_sundials
	$(CONTEXT) ./$< >$@

test_chemreac_catch.out: test_chemreac_catch
	./$< >$@

catch.hpp: catch.hpp.bz2
	bzcat $< >$@

test_chemreac_catch: test_chemreac_catch.cpp catch.hpp test_utils.o chemreac.o
	$(CXX) $(CXXFLAGS) $(FLAGS) $(INCLUDE) -o $@ $< chemreac.o test_utils.o $(LIBS)
