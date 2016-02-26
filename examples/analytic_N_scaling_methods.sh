#!/bin/bash -e
CHEMREAC_SOLVER_KWARGS="{'iterative': 0}" CHEMREAC_SOLVER=sundials python analytic_N_scaling.py -p -s sundials_direct.png
CHEMREAC_SOLVER_KWARGS="{'iterative': 1}" CHEMREAC_SOLVER=sundials python analytic_N_scaling.py -p -s sundials_iterative.png
CHEMREAC_N_JAC_DIAGS=0 CHEMREAC_SOLVER_KWARGS="{'iterative': 1}" CHEMREAC_SOLVER=sundials python analytic_N_scaling.py -p -s sundials_iterative_ndiags.png
