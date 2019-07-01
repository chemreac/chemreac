# -*- coding: utf-8 -*-
# -*- mode: cython-mode -*-

from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.unordered_map cimport unordered_map
from libcpp.string cimport string

from anyode cimport Info

cdef extern from "chemreac.hpp" namespace "chemreac":
    cdef cppclass ReactionDiffusion[T]:
        # (Private)
        T * lap_weight

        const int n, N, nr, nstencil, nsidep
        const bool logy, logt, logx, lrefl, rrefl, auto_efield
        const vector[vector[int]] stoich_active
        vector[vector[int]] stoich_inact
        const vector[vector[int]] stoich_prod
        vector[T] k
        vector[T] D
        vector[T] mobility
        vector[int] z_chg
        const vector[T] x
        vector[vector[T]] bin_k_factor
        vector[int] bin_k_factor_span
        const pair[T, T] surf_chg
        const T eps_rel
        const T faraday_const
        const T vacuum_permittivity
        vector[vector[T]] g_values
        vector[int] g_value_parents
        vector[vector[T]] fields
        vector[int] modulated_rxns
        vector[vector[T]] modulation
        vector[T] m_upper_bounds
        vector[T] m_lower_bounds
        T ilu_limit
        T m_get_dx_max_factor, m_get_dx_max_upper_limit
        T m_get_dx0_factor, m_get_dx0_max_dx
        int n_jac_diags
        bool use_log2
        bool clip_to_pos
        bool m_error_outside_bounds
        T * efield
        vector[T] gradD
        T * xc

        long nfev
        long njev
        long nprec_setup
        long nprec_solve
        long njacvec_dot
        long nprec_solve_ilu
        long nprec_solve_lu

        Info current_info
        bool autonomous_exprs, use_get_dx_max

        ReactionDiffusion(int,
                          const vector[vector[int]],
                          const vector[vector[int]],
                          vector[T],
                          int,
                          vector[T],
                          const vector[int],
                          vector[T],
                          const vector[T],
                          vector[vector[int]],
                          int,
                          bool,
                          bool,
                          bool,
                          int,
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
                          int,
                          bool,
                          bool
                          ) except +
        void zero_counters() except +
        void rhs(T, const T * const, T * const) except +
        void dense_jac_rmaj(T, const T * const, const T * const, T * const, long int) except +
        void dense_jac_cmaj(T, const T * const, const T * const, T * const, long int) except +
        void banded_jac_cmaj(T, const T * const, const T * const, T * const, long int) except +
        void compressed_jac_cmaj(T, const T * const, const T * const, T * const, long int) except +

        void per_rxn_contrib_to_fi(T, const T * const, int, T * const) except +
        int get_geom_as_int() except +
        void calc_efield(const T * const) except +

        int stencil_bi_lbound_(int) except +
        int xc_bi_map_(int) except +
