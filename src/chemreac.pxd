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
    cdef cppclass ReactionDiffusion[T]:
        # (Private)
        T * D_weight

        const uint n, N, nr, nstencil, nsidep
        const bool logy, logt, logx, lrefl, rrefl, auto_efield
        const vector[vector[uint]] stoich_active
        vector[vector[uint]] stoich_inact
        const vector[vector[uint]] stoich_prod
        vector[T] k
        vector[T] D
        vector[T] mobility
        vector[int] z_chg
        const vector[T] x
        vector[vector[T]] bin_k_factor
        vector[uint] bin_k_factor_span
        const pair[T, T] surf_chg
        const T eps_rel
        const T faraday_const
        const T vacuum_permittivity
        vector[vector[T]] g_values
        vector[int] g_value_parents
        vector[vector[T]] fields
        vector[int] modulated_rxns
        vector[vector[T]] modulation
        T ilu_limit
        uint n_jac_diags
        T * const efield
        T * xc

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
                          vector[T],
                          uint,
                          vector[T],
                          const vector[int],
                          vector[T],
                          const vector[T],
                          vector[vector[uint]],
                          int,
                          bool,
                          bool,
                          bool,
                          uint,
                          bool,
                          bool,
                          bool,
                          pair[T, T],
                          T,
                          T,
                          T,
                          vector[vector[T]],
                          vector[int],
                          vector[vector[T]],
                          vector[int],
                          vector[vector[T]],
                          T,
                          uint
                          ) except +
        void zero_counters()
        void rhs(T, const T * const, T * const)
        void dense_jac_rmaj(T, const T * const, const T * const, T * const, int)
        void dense_jac_cmaj(T, const T * const, const T * const, T * const, int)
        void banded_padded_jac_cmaj(T, const T * const, const T * const, T * const, int)
        void banded_packed_jac_cmaj(T, const T * const, const T * const, T * const, int)
        void compressed_jac_cmaj(T, const T * const, const T * const, T * const, int)

        void per_rxn_contrib_to_fi(T, const T * const, uint, T * const)
        int get_geom_as_int()
        void calc_efield(const T * const)

        uint stencil_bi_lbound_(uint)
        uint xc_bi_map_(uint)
