# -*- coding: utf-8 -*-

"""
This module is used for verification of mathematical correctness in the
implementation (used by the test suite). It uses SymPy to to derivations
symbolically. It is therefore too slow for use in integration or large systems
in general.
"""

from __future__ import (
    print_function, division, absolute_import, unicode_literals
)

from functools import reduce
from operator import add
from math import exp

import numpy as np
import sympy as sp

from .core import FLAT, CYLINDRICAL, SPHERICAL, ReactionDiffusionBase
from .util.grid import padded_centers, pxci_to_bi, stencil_pxci_lbounds


class SymRD(ReactionDiffusionBase):

    @classmethod
    def from_rd(cls, rd):
        import inspect
        return cls(*(getattr(rd, attr) for attr in
                     inspect.getargspec(cls.__init__).args[1:]))

    def __init__(self, n, stoich_reac, stoich_prod, k, N=0, D=None,
                 z_chg=None, mobility=None, x=None, stoich_actv=None,
                 bin_k_factor=None, bin_k_factor_span=None, geom=FLAT,
                 logy=False, logt=False, logx=False, nstencil=None,
                 lrefl=True, rrefl=True, auto_efield=False,
                 surf_chg=(0.0, 0.0), eps=1.0, **kwargs):
        # Save args
        self.n = n
        self.stoich_reac = stoich_reac
        self.stoich_prod = stoich_prod
        self.k = k
        self.D = D if D is not None else [0]*n
        self.z_chg = z_chg if z_chg is not None else [0]*n
        self.mobility = mobility if mobility is not None else [0]*n
        self.x = x if x is not None else [0, 1]
        self.N = len(self.x) - 1
        if N not in [None, 0]:
            assert self.N == N
        self.stoich_actv = stoich_actv or [[]*len(stoich_reac)]
        self.bin_k_factor = bin_k_factor
        self.bin_k_factor_span = bin_k_factor_span
        self.geom = geom
        self.logy = logy
        self.logt = logt
        self.logx = logx
        self.nstencil = nstencil or 3
        self.lrefl = lrefl
        self.rrefl = rrefl
        self.auto_efield = auto_efield
        self.surf_chg = surf_chg
        self.eps = eps
        if kwargs:
            raise KeyError("Don't know what to do with:", kwargs)

        # Set attributes used later
        self._t = sp.Symbol('t')
        self._nsidep = (self.nstencil-1) // 2
        self._y = sp.symbols('y:'+str(self.n*self.N))
        self._xc = padded_centers(self.x, self._nsidep)
        self._lb = stencil_pxci_lbounds(self.nstencil, self.N,
                                        self.lrefl, self.rrefl)
        self._pxci2bi = pxci_to_bi(self.nstencil, self.N)
        self._f = [0]*self.n*self.N
        self._cum_bin_k_factor_span = np.cumsum(self.bin_k_factor_span)
        self.efield = [0]*self.N

        # Reactions
        for ri, (k, sreac, sactv, sprod) in enumerate(zip(
                self.k, self.stoich_reac, self.stoich_actv,
                self.stoich_prod)):
            c_reac = map(sreac.count, range(self.n))
            c_prod = map(sprod.count, range(self.n))
            c_totl = [nprd - nrct for nrct, nprd in zip(c_reac, c_prod)]
            if sactv == []:
                sactv = sreac
            for bi in range(self.N):
                r = k * self.factor(ri, bi)
                for si in sactv:
                    r *= self.y(bi, si)
                for si in range(self.n):
                    self._f[bi*self.n+si] += c_totl[si]*r

        if self.N > 1:
            # Diffusion
            self.D_wghts = []
            self.A_wghts = []
            for bi in range(self.N):
                local_x_serie = self._xc[
                    self._lb[bi]:self._lb[bi]+self.nstencil]
                l_x_rnd = self._xc[bi+self._nsidep]
                w = sp.finite_diff_weights(2, local_x_serie, l_x_rnd)
                self.D_wghts.append(w[-1][-1])
                self.A_wghts.append(w[-2][-1])
                for wi in range(self.nstencil):
                    if self.logx:
                        if geom == FLAT:
                            self.D_wghts[bi][wi] -= w[-2][-1][wi]
                        elif geom == CYLINDRICAL:
                            self.A_wghts[bi][wi] += w[-3][-1][wi]
                        elif geom == SPHERICAL:
                            self.D_wghts[bi][wi] += w[-2][-1][wi]
                            self.A_wghts[bi][wi] += 2*w[-3][-1][wi]
                        self.D_wghts[bi][wi] *= exp(-2*l_x_rnd)
                        self.A_wghts[bi][wi] *= exp(-l_x_rnd)
                    else:
                        if geom == CYLINDRICAL:
                            self.D_wghts[bi][wi] += w[-2][-1][wi]/l_x_rnd
                            self.A_wghts[bi][wi] += w[-3][-1][wi]/l_x_rnd
                        elif geom == SPHERICAL:
                            self.D_wghts[bi][wi] += 2*w[-2][-1][wi]/l_x_rnd
                            self.A_wghts[bi][wi] += 2*w[-3][-1][wi]/l_x_rnd

            for bi, (dw, aw) in enumerate(zip(self.D_wghts, self.A_wghts)):
                for si in range(self.n):
                    d_terms = [dw[k]*self.y(
                        self._pxci2bi[self._lb[bi]+k], si
                    ) for k in range(self.nstencil)]
                    self._f[bi*self.n + si] += self.D[si]*reduce(
                        add, d_terms)
                    a_terms = [aw[k]*self.y(
                        self._pxci2bi[self._lb[bi]+k], si
                    ) for k in range(self.nstencil)]
                    self._f[bi*self.n + si] += (
                        self.mobility[si]*self.efield[bi]*reduce(
                            add, a_terms))

        if self.logy or self.logt:
            for bi in range(self.N):
                for si in range(self.n):
                    if self.logy:
                        self._f[bi*self.n+si] /= self.y(bi, si)
                    if self.logt:
                        self._f[bi*self.n+si] *= sp.exp(self._t)

    def factor(self, ri, bi):
        for i, ub in enumerate(self._cum_bin_k_factor_span):
            if ri < ub:
                idx = i
                break
        else:
            return 1.0
        return self.bin_k_factor[bi][idx]

    def y(self, bi, si):
        if self.logy:
            return sp.exp(self._y[bi*self.n+si])
        else:
            return self._y[bi*self.n+si]

    def f(self, t, y, fout):
        subsd = dict(zip(self._y, y))
        subsd[self._t] = t
        fout[:] = [expr.subs(subsd) for expr in self._f]

    @property
    def jacobian(self):
        try:
            return self._jacobian
        except AttributeError:
            fmat = sp.Matrix(1, self.n*self.N, lambda q, i: self._f[i])
            self._jacobian = fmat.jacobian(self._y)
            return self._jacobian

    def dense_jac_rmaj(self, t, y, Jout):
        subsd = dict(zip(self._y, y))
        subsd[self._t] = t
        Jout[:, :] = [[expr.subs(subsd) for expr in row]
                      for row in self.jacobian.tolist()]
