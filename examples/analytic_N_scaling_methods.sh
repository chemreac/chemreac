#!/bin/bash -ex
# ./analytic_N_scaling_methods.sh --ylims '0.03125,33554432;10,1e5;0,1e4;1e-4,1e2' --nfit '7,5' --nNs 4
DEST=./analytic_N_scaling_output
for METH in adams bdf; do
    echo "DIRECT NDIAGS=1 $METH"
    CHEMREAC_N_JAC_DIAGS=1 CHEMREAC_INTEGRATION_KWARGS="{'linear_solver': 'banded', 'method': '$METH'}" CHEMREAC_SOLVER=sundials ${PYTHON:-python} analytic_N_scaling.py --plot --savefig $DEST/${METH}_sundials_direct.png $@
    echo "DIRECT NDIAGS=0 $METH"
    CHEMREAC_N_JAC_DIAGS=0 CHEMREAC_INTEGRATION_KWARGS="{'linear_solver': 'banded'}" CHEMREAC_SOLVER=sundials ${PYTHON:-python} analytic_N_scaling.py --plot --savefig $DEST/${METH}_sundials_direct_ndiags.png $@
    echo "ITERATIVE NDIAGS=1 $METH"
    CHEMREAC_N_JAC_DIAGS=1 CHEMREAC_INTEGRATION_KWARGS="{'linear_solver': 'gmres', 'method': '$METH'}" CHEMREAC_SOLVER=sundials ${PYTHON:-python} analytic_N_scaling.py --plot --savefig $DEST/${METH}_sundials_iterative.png $@
    echo "ITERATIVE NDIAGS=0 $METH"
    CHEMREAC_N_JAC_DIAGS=0 CHEMREAC_INTEGRATION_KWARGS="{'linear_solver': 'gmres', 'method': '$METH'}" CHEMREAC_SOLVER=sundials python analytic_N_scaling.py --plot --savefig $DEST/${METH}_sundials_iterative_ndiags.png $@
done
