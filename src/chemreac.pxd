# -*- coding: utf-8 -*-
# -*- mode: cython-mode -*-

from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.unordered_map cimport unordered_map
from libcpp.string cimport string

cdef extern from *:
    ctypedef unsigned int uint

cdef extern from "chemreac.hpp" namespace "chemreac":
    cdef cppclass ReactionDiffusion:
        # (Private)
        double * D_weight

        const uint n, N, nr, nstencil, nsidep
        const bool logy, logt, logx, lrefl, rrefl, auto_efield
        const vector[vector[uint]] stoich_active
        vector[vector[uint]] stoich_inact
        const vector[vector[uint]] stoich_prod
        vector[double] k
        vector[double] D
        vector[double] mobility
        vector[int] z_chg
        const vector[double] x
        vector[vector[double]] bin_k_factor
        vector[uint] bin_k_factor_span
        const pair[double, double] surf_chg
        const double eps_rel
        vector[vector[double]] g_values
        vector[int] g_value_parents
        vector[vector[double]] fields
        vector[int] modulated_rxns
        vector[vector[double]] modulation
        double ilu_limit
        uint n_jac_diags
        double * const efield
        double * xc

        long nfev
        long njev
        long nprec_setup
        long nprec_solve
        long njacvec_dot
        long nprec_solve_ilu
        long nprec_solve_lu

        unordered_map[string, int] last_integration_info

        ReactionDiffusion(uint,
                          const vector[vector[uint]],
                          const vector[vector[uint]],
                          vector[double],
                          uint,
                          vector[double],
                          const vector[int],
                          vector[double],
                          const vector[double],
                          vector[vector[uint]],
                          int,
                          bool,
                          bool,
                          bool,
                          uint,
                          bool,
                          bool,
                          bool,
                          pair[double, double],
                          double,
                          double,
                          double,
                          vector[vector[double]],
                          vector[int],
                          vector[vector[double]],
                          vector[int],
                          vector[vector[double]],
                          double,
                          uint
                          ) except +
        void zero_counters()
        void f(double, const double * const, double * const)
        void dense_jac_rmaj(double, const double * const, const double * const, double * const, int)
        void dense_jac_cmaj(double, const double * const, const double * const, double * const, int)
        void banded_padded_jac_cmaj(double, const double * const, const double * const, double * const, int)
        void banded_packed_jac_cmaj(double, const double * const, const double * const, double * const, int)
        void compressed_jac_cmaj(double, const double * const, const double * const, double * const, int)

        void per_rxn_contrib_to_fi(double, const double * const, uint, double * const)
        int get_geom_as_int()
        void calc_efield(const double * const)

        uint stencil_bi_lbound_(uint)
        uint xc_bi_map_(uint)
