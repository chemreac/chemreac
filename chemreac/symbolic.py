# -*- coding: utf-8 -*-

from __future__ import print_function, division, absolute_import, unicode_literals

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

    def __init__(self, n, stoich_reac, stoich_prod, k, N=0, D=None, x=None,
                stoich_actv=None, bin_k_factor=None, bin_k_factor_span=None,
                geom=FLAT, logy=False, logt=False, nstencil=None, lrefl=True,
                rrefl=True, **kwargs):
        # Save args
        self.n = n
        self.stoich_reac = stoich_reac
        self.stoich_prod = stoich_prod
        self.k = k
        self.D = D if D != None else [0]*n
        self.x = x if x != None else [0, 1]
        self.N = N if N != None else len(self.x) - 1
        self.stoich_actv = stoich_actv or [[]*len(stoich_reac)]
        self.bin_k_factor = bin_k_factor
        self.bin_k_factor_span = bin_k_factor_span
        self.geom = geom
        self.logy = logy
        self.logt = logt
        self.nstencil = nstencil or 3
        self.lrefl = lrefl
        self.rrefl = rrefl
        if kwargs: raise KeyError("Don't know what to do with:", kwargs)

        # Set attributes used later
        self._t = sp.Symbol('t')
        self._nsidep = (self.nstencil-1) // 2
        self._y = sp.symbols('y:'+str(self.n*self.N))
        self._xc = padded_centers(self.x, self._nsidep)
        self._lb = stencil_pxci_lbounds(self.nstencil, self.N, self.lrefl, self.rrefl)
        self._pxci2bi = pxci_to_bi(self.nstencil, self.N)
        self._f = [0]*self.n*self.N
        self._cum_bin_k_factor_span = np.cumsum(self.bin_k_factor_span)

        # Reactions
        for ri, (k, sreac, sactv, sprod) in enumerate(zip(
                self.k, self.stoich_reac, self.stoich_actv, self.stoich_prod)):
            c_reac = [sreac.count(i) for i in range(self.n)]
            c_prod = [sprod.count(i) for i in range(self.n)]
            c_totl = [nprod - nreac for nreac, nprod in zip(c_reac, c_prod)]
            if sactv == []:
                sactv = sreac
            for bi in range(self.N):
                r = k * self.factor(ri, bi)
                for si in sactv:
                    r *= self.y(bi, si)
                for si in range(self.n):
                    self._f[bi*self.n+si] += c_totl[si]*r

        print(geom)
        if self.N > 1:
            # Diffusion
            self.D_weights = []
            for bi in range(self.N):
                local_x_serie = self._xc[self._lb[bi]:self._lb[bi]+self.nstencil]
                local_x_around = self._xc[bi+self._nsidep]
                w = sp.finite_diff_weights(2, local_x_serie, local_x_around)
                self.D_weights.append(w[-1][-1])
                for wi in range(self.nstencil):
                    if geom == CYLINDRICAL:
                        self.D_weights[bi][wi] += w[-2][-1][wi]/local_x_around
                    if geom == SPHERICAL:
                        self.D_weights[bi][wi] += 2*w[-2][-1][wi]/local_x_around
            for bi, w in enumerate(self.D_weights):
                for si in range(self.n):
                    fd_terms = [w[k]*self.y(self._pxci2bi[self._lb[bi]+k], si)
                                for k in range(self.nstencil)]
                    self._f[bi*self.n + si] += self.D[si]*reduce(add, fd_terms)

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

    def dense_jac_rmaj(self, t, y, Jout):
        subsd = dict(zip(self._y, y))
        subsd[self._t] = t
        fmat = sp.Matrix(1, self.n*self.N, lambda q, i: self._f[i])
        Jout[:,:] = [[expr.subs(subsd) for expr in row]
                     for row in fmat.jacobian(self._y).tolist()]
