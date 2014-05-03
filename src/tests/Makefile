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
  CXXFLAGS=-std=c++0x -DDEBUG
  OPENMPLIBS=-liomp5
  OPENMPFLAG=-openmp
  C99=-std=c99
else
ifeq ($(COMPILER_VENDOR),cray-pgi)
  FC=ftn
  CXX=CC
  CC=cc
  FORTCPPLINKER=CC
  FORTCPPLIBS=-pgf90libs
  FORTCPPFLAGS=
  FORTFLAGS=
  PERFFLAGS=-O2
  FLAGS=-DCRAY
  CXXFLAGS= -DDEBUG

  C99=-c99
  OPENMPLIBS=
  OPENMPFLAG=-mp
else
  PROFFLAGS=-pg
ifeq ($(COMPILER_VENDOR),cray-gnu)
  FC=ftn
  CXX=CC
  CC=cc
  FORTCPPLINKER=CC
  FORTCPPLIBS=-lgfortranbegin -lgfortran
  FORTCPPFLAGS=
  FORTFLAGS=-std=f2008 
  PERFFLAGS=-funroll-loops -O3 -march=native -ffast-math
  FLAGS=-Wall -Wextra -pedantic -DCRAY
  CXXFLAGS=-std=c++0x -DDEBUG
  C99=-std=c99
  OPENMPLIBS=#-lgomp
  OPENMPFLAG=-fopenmp
else
  FC=gfortran
  CXX=g++ 
  CC=gcc
  FORTCPPLINKER=g++
  FORTCPPLIBS=-lgfortranbegin -lgfortran
  FORTCPPFLAGS=
  FORTFLAGS=-std=f2008 
  PERFFLAGS=-funroll-loops -O3 -march=native -ffast-math
  FLAGS=-Wall -Wextra -pedantic
  OPENMPLIBS=-lgomp
  OPENMPFLAG=-fopenmp
  CXXFLAGS=-std=c++0x -g -DDEBUG # $(OPENMPFLAG) $(PERFFLAGS)
  LIBS=-lrt # $(OPENMPLIBS)
  C99=-std=c99
endif
endif
endif

ifeq ($(OPTIMIZE),1)
  FLAGS += $(PERFFLAGS)
  FORTFLAGS += $(PERFFLAGS)
else
  FLAGS +=-g -DDEBUG 
  FORTFLAGS +=-g #-DDEBUG
  ifeq ($(COMPILER_VENDOR),intel) # 2013 SP 1
    FLAGS +=-debug all
    FORTFLAGS +=-debug all
  endif
endif
ifeq ($(PROFILING),1)
  FLAGS += $(PROFFLAGS) -DPROFILING 
  FORTFLAGS += $(PROFFLAGS) -DPROFILING 
endif

FLAGS += $(EXTRA_FLAGS)

#note if missing libraries do, e.g.:
#env LD_LIBRARY_PATH=/usr/lib/gcc/x86_64-linux-gnu/4.6/:$LD_LIBRARY_PATH make

%.o: %.cpp
	$(CXX) -fPIC $(CXXFLAGS) $(FLAGS) -c -o $@ $<


all: test_chemreac.out

.PHONY: clean all

GENERATED=chemreac.cpp test_chemreac test_chemreac.out fornberg.o

clean:
	rm $(GENERATED)

chemreac.cpp: ../chemreac_template.cpp
	python enmako.py -g '{"USE_OPENMP": False}' -o $@ $<

fornberg.o: ../finitediff/finitediff/fornberg.f90
	$(FC) $(FORTFLAGS) -c -o $@ $^

c_fornberg.o: ../finitediff/finitediff/c_fornberg.f90
	$(FC) $(FORTFLAGS) -c -o $@ $^


test_chemreac: chemreac.cpp test_chemreac.cpp fornberg.o c_fornberg.o
	$(CXX) $(PROFFLAGS) $(CXXFLAGS) -I.. -I../finitediff/finitediff -Drestrict=__restrict__\
	 -o $@ $^ $(LIBS) $(FORTCPPLIBS) # -lrt only needed if we're not using OpenMP

test_chemreac.out: test_chemreac
	./$< >$@

profile_test_chemreac.out: test_chemreac
	./$<
	gprof ./$< > $@
