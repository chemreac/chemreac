#!/bin/bash -ex
# ./analytic_N_scaling_methods.sh --ylims '0.03125,33554432;10,1e5;0,1e4;1e-4,1e2' --nfit '7,5' --nNs 4
DEST=./analytic_N_scaling_output
for METH in adams bdf; do
    echo "DIRECT NDIAGS=1 $METH"
    CHEMREAC_N_JAC_DIAGS=1 CHEMREAC_SOLVER_KWARGS="{'iterative': 0}" CHEMREAC_SOLVER=sundials ${PYTHON:-python} analytic_N_scaling.py --meth $METH --plot --savefig $DEST/sundials_direct_$METH.png $@
    echo "DIRECT NDIAGS=0 $METH"
    CHEMREAC_N_JAC_DIAGS=0 CHEMREAC_SOLVER_KWARGS="{'iterative': 0}" CHEMREAC_SOLVER=sundials ${PYTHON:-python} analytic_N_scaling.py --meth $METH --plot --savefig $DEST/sundials_direct_ndiags_$METH.png $@
    echo "ITERATIVE NDIAGS=1 $METH"
    CHEMREAC_N_JAC_DIAGS=1 CHEMREAC_SOLVER_KWARGS="{'iterative': 1, 'method': '$METH'}" CHEMREAC_SOLVER=sundials ${PYTHON:-python} analytic_N_scaling.py --plot --savefig $DEST/sundials_iterative_$METH.png $@
    echo "ITERATIVE NDIAGS=0 $METH"
    CHEMREAC_N_JAC_DIAGS=0 CHEMREAC_SOLVER_KWARGS="{'iterative': 1, 'method': '$METH'}" CHEMREAC_SOLVER=sundials python analytic_N_scaling.py --plot --savefig $DEST/sundials_iterative_ndiags_$METH.png $@
done
