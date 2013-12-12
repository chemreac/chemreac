## When using intel compiler: source e.g. /opt/intel/mkl/bin/mklvars.sh to get $MKLROOT

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
  CXXFLAGS=-std=c++0x
  OPENMPLIBS=-liomp5
  OPENMPFLAG=-openmp
  LAPACKLIBS=$(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a -lpthread -lm 
  LAPACKFLAGS=-openmp -I$(MKLROOT)/include/intel64/lp64 -mkl=parallel
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
  CXXFLAGS=
  LAPACKLIBS=
  LAPACKFLAGS=
  C99=-c99
  OPENMPLIBS=
  OPENMPFLAG=-mp
else
  PROFFLAGS=-pg
ifeq ($(COMPILER_VENDOR),cray-gnu)
#-lpgf90 -lpgf90_rpm1 -lpgf902 -lpgf90rtl -lpgftnrtl
  FC=ftn
  CXX=CC
  CC=cc
  FORTCPPLINKER=CC
  FORTCPPLIBS=-lgfortranbegin -lgfortran
  FORTCPPFLAGS=
  FORTFLAGS=-std=f2008 
  PERFFLAGS=-funroll-loops -O3 -march=native -ffast-math
  FLAGS=-Wall -Wextra -pedantic -DCRAY
  CXXFLAGS=-std=c++0x
  LAPACKLIBS=
  LAPACKFLAGS=
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
  CXXFLAGS=-std=c++0x
  LAPACKLIBS=-llapack -lblas
  LAPACKFLAGS=
  C99=-std=c99
  OPENMPLIBS=-lgomp
  OPENMPFLAG=-fopenmp
endif
endif
endif
EXTRA_FLAGS=

ifeq ($(OPTIMIZE),1)
  FLAGS += $(PERFFLAGS)
  FORTFLAGS += $(PERFFLAGS)
else
  FLAGS +=-g #-DDEBUG 
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

%.o: %.f90
	$(FC) -fPIC $(FORTFLAGS) $(FLAGS) -c -o $@ $<

all: test_cpp_chem

.PHONY: clean all

clean:
	rm $(GENERATED)
	rm $(OBJS)

cpp_chem.cpp: cpp_chem_template.cpp
	python enmako.py -g '{}' -o $@ $<

test_cpp_chem: cpp_chem.cpp test_cpp_chem.cpp
	$(CXX) $(PERFFLAGS) $(PROFFLAGS) $(CXXFLAGS) ${LAPACKFLAGS} -o $@ $^ -lrt $(LAPACKLIBS) $(OPENMPLIBS)

profile_test_cpp_chem.out: test_cpp_chem
	./$<
	gprof ./$< > $@
