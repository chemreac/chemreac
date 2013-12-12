hpcprj
======

Chemical Reaction-Diffusion ODE systems discretized in
one dimension. The Jacobian is used by implicit methods (which is usually needed
since the systems are stiff)

Setup
=====
``` python setup.py build ``` is enough on most *NIX machines.

Parallelized evaluation of f and jac
====================================
```
env OMP_NUM_THREADS=1 python -m cProfile -o profile.out radiolysis.py -t 60 -N 300 -i cpp && runsnake profile.out
```


Installation
============
To install using Intel's compilers and Intel MKL, you need to set environment
variables. On my system (MKL 11.1) I would do:

source /opt/intel/mkl/bin/mklvars.sh intel64 lp64; COMPILER_VENDOR=intel python setup.py

for i in {1,2,3,4}; do env OMP_NUM_THREADS=$i python -m cProfile -o profile_${i}.out radiolysis.py -i cpp -N 3000; done

