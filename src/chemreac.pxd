# -*- coding: utf-8 -*-
# -*- mode: cython-mode -*-

from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.utility cimport pair

cdef extern from *:
    ctypedef unsigned int uint

cdef extern from "chemreac.hpp" namespace "chemreac":
    cdef cppclass ReactionDiffusion:
        # (Private)
        double * D_weight

        const uint n, N, nr, nstencil
        const bool logy, logt, logx, lrefl, rrefl, auto_efield
        const vector[vector[uint]] stoich_reac
        vector[vector[uint]] stoich_actv
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
        double * const efield
        double * xc
        long neval_f
        long neval_j
        long nprec_setup
        long nprec_solve
        long njacvec_dot

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
                          vector[vector[double]]
                          ) except +
        void f(double, const double * const, double * const)
        void dense_jac_rmaj(double, const double * const, const double * const, double * const, int)
        void dense_jac_cmaj(double, const double * const, const double * const, double * const, int)
        void banded_padded_jac_cmaj(double, const double * const, const double * const, double * const, int)
        void banded_packed_jac_cmaj(double, const double * const, const double * const, double * const, int)
        void compressed_jac_cmaj(double, const double * const, const double * const, double * const, int)

        void per_rxn_contrib_to_fi(double, const double * const, uint, double * const)
        int get_geom_as_int()
        void calc_efield(const double * const)

        uint _stencil_bi_lbound(uint)
        uint _xc_bi_map(uint)
