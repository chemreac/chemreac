# -*- coding: utf-8 -*-
# -*- mode: cython-mode -*-

from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.utility cimport pair

cdef extern from *:
    ctypedef unsigned int uint

cdef extern from "chemreac.h" namespace "chemreac":
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
        const double eps
        double * xc
        double * const efield
        uint neval_f
        uint neval_j


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
                          vector[vector[double]],
                          vector[uint],
                          int,
                          bool,
                          bool,
                          bool,
                          uint,
                          bool,
                          bool,
                          bool,
                          pair[double, double],
                          double) except +
        void f(double, const double * const, double * const)
        void dense_jac_rmaj(double, const double * const, double * const, int)
        void dense_jac_cmaj(double, const double * const, double * const, int)
        void banded_padded_jac_cmaj(double, const double * const, double * const, int)
        void banded_packed_jac_cmaj(double, const double * const, double * const, int)

        void per_rxn_contrib_to_fi(double, const double * const, uint, double * const)
        int get_geom_as_int()
        void calc_efield(const double * const)

        uint _stencil_bi_lbound(uint)
        uint _xc_bi_map(uint)
