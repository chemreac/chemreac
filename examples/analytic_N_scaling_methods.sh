#!/bin/bash -e
# ./analytic_N_scaling_methods.sh --ylims '0.03125,33554432;10,1e5;0,1e4;1e-4,1e2' --nfit '7,5' --nNs 4
for METH in adams bdf; do
    CHEMREAC_SOLVER_KWARGS="{'iterative': 0, 'method': '$METH'}" CHEMREAC_SOLVER=sundials python analytic_N_scaling.py -p -s sundials_direct_$METH.png $@
    CHEMREAC_SOLVER_KWARGS="{'iterative': 1, 'method': '$METH'}" CHEMREAC_SOLVER=sundials python analytic_N_scaling.py -p -s sundials_iterative_$METH.png $@
    CHEMREAC_N_JAC_DIAGS=0 CHEMREAC_SOLVER_KWARGS="{'iterative': 1, 'method': '$METH'}" CHEMREAC_SOLVER=sundials python analytic_N_scaling.py -p -s sundials_iterative_ndiags_$METH.png $@
done
