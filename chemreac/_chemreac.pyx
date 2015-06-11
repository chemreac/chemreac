# -*- coding: utf-8 -*-
# distutils: language = c++

import cython

import numpy as np
cimport numpy as cnp

from chemreac cimport ReactionDiffusion
from cvodes_wrapper cimport simple_integrate

from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.utility cimport pair

cdef extern from *:
    ctypedef unsigned int uint

DEF FLAT=0
DEF CYLINDRICAL=1
DEF SPHERICAL=2


cdef class ArrayWrapper(object):
    cdef public dict __array_interface__

    def __init__(self, **kwargs):
        self.__array_interface__ = kwargs


cdef fromaddress(address, shape, dtype=np.float64, strides=None, ro=True):
    dtype = np.dtype(dtype)
    return np.asarray(ArrayWrapper(
        data=(address, ro),
        typestr=dtype.str,
        descr=dtype.descr,
        shape=shape,
        strides=strides,
        version=3,
    ))

def diag_data_len(N, n, ndiag):
    return n*(N*ndiag - ((ndiag+1)*(ndiag+1) - (ndiag+1))//2)


cdef class CppReactionDiffusion:
    """
    Wrapper around C++ class ReactionDiffusion,
    """
    cdef ReactionDiffusion *thisptr
    cdef public vector[double] k_err, D_err
    cdef public list names, tex_names

    def __cinit__(self,
                  uint n,
                  vector[vector[uint]] stoich_active,
                  vector[vector[uint]] stoich_prod,
                  vector[double] k,
                  uint N,
                  vector[double] D,
                  vector[int] z_chg,
                  vector[double] mobility,
                  vector[double] x,
                  vector[vector[uint]] stoich_inactv,
                  int geom,
                  bint logy,
                  bint logt,
                  bint logx,
                  vector[vector[double]] g_values,
                  vector[int] g_value_parents,
                  vector[vector[double]] fields,
                  vector[int] modulated_rxns,
                  vector[vector[double]] modulation,
                  uint nstencil=3,
                  bint lrefl=True,
                  bint rrefl=True,
                  bint auto_efield=False,
                  pair[double, double] surf_chg=(0, 0),
                  double eps_rel=1.0,
                  double faraday_const=9.64853399e4,
                  double vacuum_permittivity=8.854187817e-12,
              ):
        cdef size_t i
        self.thisptr = new ReactionDiffusion(
            n, stoich_active, stoich_prod, k, N,
            D, z_chg, mobility, x, stoich_inactv, geom,
            logy, logt, logx, nstencil,
            lrefl, rrefl, auto_efield, surf_chg, eps_rel, faraday_const,
            vacuum_permittivity, g_values, g_value_parents, fields,
            modulated_rxns, modulation)

    def __dealloc__(self):
        del self.thisptr

    def f(self, double t, cnp.ndarray[cnp.float64_t, ndim=1] y,
          cnp.ndarray[cnp.float64_t, ndim=1] fout):
        assert y.size == fout.size
        assert y.size >= self.n
        self.thisptr.f(t, &y[0], &fout[0])

    def dense_jac_rmaj(self, double t, cnp.ndarray[cnp.float64_t, ndim=1] y,
                       cnp.ndarray[cnp.float64_t, ndim=2, mode="c"] Jout):
        assert y.size >= self.n*self.N
        assert Jout.shape[0] >= self.n*self.N
        assert Jout.shape[1] >= self.n*self.N
        self.thisptr.dense_jac_rmaj(
            t, &y[0], NULL, &Jout[0, 0], Jout.shape[1])

    def dense_jac_cmaj(self, double t, cnp.ndarray[cnp.float64_t, ndim=1] y,
                       cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran"] Jout):
        assert y.size >= self.n*self.N
        assert Jout.shape[0] >= self.n*self.N
        assert Jout.shape[1] >= self.n*self.N
        self.thisptr.dense_jac_cmaj(
            t, &y[0], NULL, &Jout[0, 0], Jout.shape[0])

    def banded_padded_jac_cmaj(self, double t, cnp.ndarray[cnp.float64_t, ndim=1] y,
                       cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran"] Jout):
        assert y.size >= self.n*self.N
        assert Jout.shape[0] >= self.n*3+1
        assert Jout.shape[1] >= self.n*self.N
        self.thisptr.banded_padded_jac_cmaj(
            t, &y[0], NULL, &Jout[0, 0], Jout.shape[0])

    def banded_packed_jac_cmaj(self, double t, cnp.ndarray[cnp.float64_t, ndim=1] y,
                       cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran"] Jout):
        assert y.size >= self.n*self.N
        assert Jout.shape[0] >= self.n*2+1
        assert Jout.shape[1] >= self.n*self.N
        self.thisptr.banded_packed_jac_cmaj(
            t, &y[0], NULL, &Jout[0, 0], Jout.shape[0])

    def compressed_jac_cmaj(self, double t, cnp.ndarray[cnp.float64_t, ndim=1] y,
                            cnp.ndarray[cnp.float64_t, ndim=1] Jout):
        assert y.size >= self.n*self.N
        assert Jout.size >= self.n*self.n*self.N + 2*diag_data_len(
            self.N, self.n, (self.nstencil-1)//2)
        self.thisptr.compressed_jac_cmaj(
            t, &y[0], NULL, <double *>Jout.data, self.n)

    def calc_efield(self, cnp.ndarray[cnp.float64_t, ndim=1] linC):
        self.thisptr.calc_efield(&linC[0])
        return self.efield  # convenience

    def integrated_conc(self, linC):
        """
        Integrates the concentration over the volume of the system.
        Pass linear concentration "linC"
        """
        if linC.shape != (self.N,):
            raise ValueError("linC must be of length N")
        if self.geom == FLAT:
            return np.sum(np.diff(self.lin_x)*linC)
        elif self.geom == CYLINDRICAL:
            return np.sum(np.pi*np.diff(self.lin_x**2)*linC)
        elif self.geom == SPHERICAL:
            return np.sum(4*np.pi/3*np.diff(self.lin_x**3)*linC)

    property lin_x:
        def __get__(self):
            if self.logx:
                return np.exp(self.x)
            else:
                return self.x

    property n:
        def __get__(self):
            return self.thisptr.n

    property N:
        def __get__(self):
            return self.thisptr.N

    property nr:
        def __get__(self):
            return self.thisptr.nr

    property geom:
        def __get__(self):
            return self.thisptr.get_geom_as_int()

    property stoich_active:
        def __get__(self):
            return self.thisptr.stoich_active

    property stoich_prod:
        def __get__(self):
            return self.thisptr.stoich_prod

    property stoich_inactv:
        def __get__(self):
            return self.thisptr.stoich_inactv

    property _k:
        def __get__(self):
            return np.asarray(self.thisptr.k)
        def __set__(self, vector[double] k):
            assert len(k) == self.nr
            self.thisptr.k = k

    property _D:
        def __get__(self):
            return np.asarray(self.thisptr.D)

        def __set__(self, vector[double] D):
            cdef size_t i
            assert len(D) == self.n
            self.thisptr.D = D

    property z_chg:
        def __get__(self):
            return np.asarray(self.thisptr.z_chg, dtype=np.int32)

        def __set__(self, vector[int] z_chg):
            assert len(z_chg) == self.n
            self.thisptr.z_chg = z_chg

    property _mobility:
        def __get__(self):
            return np.asarray(self.thisptr.mobility)

        def __set__(self, vector[double] mobility):
            cdef size_t i
            assert len(mobility) == self.n
            self.thisptr.mobility = mobility

    property x:
        def __get__(self):
            return np.asarray(self.thisptr.x)

    property logy:
        def __get__(self):
            return self.thisptr.logy

    property logt:
        def __get__(self):
            return self.thisptr.logt

    property logx:
        def __get__(self):
            return self.thisptr.logx

    property nstencil:
        def __get__(self):
            return self.thisptr.nstencil

    property lrefl:
        def __get__(self):
            return self.thisptr.lrefl

    property rrefl:
        def __get__(self):
            return self.thisptr.rrefl

    property auto_efield:
        def __get__(self):
            return self.thisptr.auto_efield

    property surf_chg:
        def __get__(self):
            return self.thisptr.surf_chg

    property eps_rel:
        def __get__(self):
            return self.thisptr.eps_rel

    property _g_values:
        def __get__(self):
            return self.thisptr.g_values
        def __set__(self, vector[vector[double]] g_values):
            self.thisptr.g_values = g_values

    property g_value_parents:
        def __get__(self):
            return self.thisptr.g_value_parents
        def __set__(self, vector[int] g_value_parents):
            self.thisptr.g_value_parents = g_value_parents

    property _fields:
        def __get__(self):
            return self.thisptr.fields
        def __set__(self, vector[vector[double]] fields):
            assert len(fields) == len(self.g_values)
            self.thisptr.fields = fields

    property modulated_rxns:
        def __get__(self):
            return self.thisptr.modulated_rxns
        def __set__(self, vector[int] modulated_rxns):
            self.thisptr.modulated_rxns = modulated_rxns

    property modulation:
        def __get__(self):
            return self.thisptr.modulation
        def __set__(self, vector[vector[double]] modulation):
            self.thisptr.modulation = modulation

    property neval_f:
        def __get__(self):
            return self.thisptr.neval_f
        def __set__(self, long n):
            self.thisptr.neval_f = n

    property neval_j:
        def __get__(self):
            return self.thisptr.neval_j
        def __set__(self, long n):
            self.thisptr.neval_j = n

    property nprec_setup:
        def __get__(self):
            return self.thisptr.nprec_setup
        def __set__(self, long n):
            self.thisptr.nprec_setup = n

    property nprec_solve:
        def __get__(self):
            return self.thisptr.nprec_solve
        def __set__(self, long n):
            self.thisptr.nprec_solve = n

    property njacvec_dot:
        def __get__(self):
            return self.thisptr.njacvec_dot
        def __set__(self, long n):
            self.thisptr.njacvec_dot = n

    def zero_out_counters(self):
        self.neval_f = 0
        self.neval_j = 0
        self.nprec_setup = 0
        self.nprec_solve = 0
        self.njacvec_dot = 0

    # Extra convenience
    def per_rxn_contrib_to_fi(self, double t, cnp.ndarray[cnp.float64_t, ndim=1] y,
                              int si, cnp.ndarray[cnp.float64_t, ndim=1] out):
        """
        Decomposes the rate of change of a species concentration into per
        reaction contributions.

        Parameters
        ----------
        t: float
            time
        y: 1-dimensional numpy array
            dependent variables
        si: int
            specie index to analyse for
        out: 1-dimensional numpy array (slice)
            output argument, of length ``nr`` (per reaction contribution)

        Returns
        -------
        None, see ``out`` parameter
        """
        self.thisptr.per_rxn_contrib_to_fi(t, &y[0], si, &out[0])

    property xcenters:
        def __get__(self):
            return fromaddress(<long>(
                &self.thisptr.xc[(self.thisptr.nstencil-1)//2]), (self.N,))

    # For debugging
    property _xc:
        def __get__(self):
            return fromaddress(<long>self.thisptr.xc,
                               (self.N + self.thisptr.nstencil - 1,))

    property D_weight:  # (Private)
        def __get__(self):
            return fromaddress(<long>self.thisptr.D_weight,
                               (self.thisptr.N*self.thisptr.nstencil,))

    def _stencil_bi_lbound(self, uint bi):
        return self.thisptr._stencil_bi_lbound(bi)

    def _xc_bi_map(self, uint xci):
        return self.thisptr._xc_bi_map(xci)

    property efield:
        def __get__(self):
            return fromaddress(<long>self.thisptr.efield, (self.thisptr.N,))
        def __set__(self, double[:] efield):
            cdef int i
            assert efield.size == self.thisptr.N
            for i in range(self.thisptr.N):
                self.thisptr.efield[i] = efield[i]

# sundials wrapper:

def sundials_integrate(
        CppReactionDiffusion rd, cnp.ndarray[cnp.float64_t, ndim=1] y0,
        cnp.ndarray[cnp.float64_t, ndim=1] tout,
        vector[double] atol, double rtol, basestring method, bool with_jacobian=True,
        int iterative=0):
    cdef cnp.ndarray[cnp.float64_t, ndim=1] yout = np.empty(tout.size*rd.n*rd.N)
    assert y0.size == rd.n*rd.N
    simple_integrate[double, ReactionDiffusion](
        rd.thisptr, atol, rtol, {'adams': 1, 'bdf': 2}[method.lower()],
        &y0[0], tout.size, &tout[0], &yout[0], with_jacobian, iterative)
    return yout.reshape((tout.size, rd.N, rd.n))


# Below is an implementation of the classic Runge Kutta 4th order stepper with fixed step size
# it is only useful for debugging purposes (fixed step size isn't for production runs)

cdef void _add_2_vecs(int n,  double * v1, double * v2,
                      double f1, double f2, double * out):
    cdef int i
    for i in range(n):
        out[i] = f1*v1[i] + f2*v2[i]


cdef void _add_5_vecs(int n, double * v1, double * v2, double * v3,
                      double * v4, double * v5, double f1, double f2,
                      double f3, double f4, double f5, double * out):
    cdef int i
    for i in range(n):
        out[i] = f1*v1[i] + f2*v2[i] + f3*v3[i] + f4*v4[i] + f5*v5[i]


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void _rk4(ReactionDiffusion * rd,
               cnp.ndarray[cnp.float64_t, ndim=1, mode='c'] y0,
               cnp.ndarray[cnp.float64_t, ndim=1, mode='c'] tout,
               cnp.ndarray[cnp.float64_t, ndim=2, mode='c'] y0out,
               cnp.ndarray[cnp.float64_t, ndim=2, mode='c'] y1out):
    # Runge-Kutta 4th order stepper
    # see: http://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods
    cdef int i
    cdef double t, h
    cdef int ny = y0.size
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode='c'] tmp = np.empty(ny, dtype=np.float64)
    cdef double *k1
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode='c'] k2 = np.empty(ny, dtype=np.float64)
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode='c'] k3 = np.empty(ny, dtype=np.float64)
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode='c'] k4 = np.empty(ny, dtype=np.float64)
    y0out[0, :] = y0[:]
    for i in range(1, tout.size):
        t = tout[i]
        h = t - tout[i-1]
        k1 = &y1out[i-1, 0]
        rd.f(t, &y0out[i-1, 0], &k1[0])
        _add_2_vecs(ny, &y0out[i-1, 0], &k1[0], 1.0, h/2, &tmp[0])
        rd.f(t + h/2, &tmp[0], &k2[0])
        _add_2_vecs(ny, &y0out[i-1, 0], &k2[0], 1.0, h/2, &tmp[0])
        rd.f(t + h/2, &tmp[0], &k3[0])
        _add_2_vecs(ny, &y0out[i-1, 0], &k3[0], 1.0, h/2, &tmp[0])
        rd.f(t + h, &tmp[0], &k4[0])
        _add_5_vecs(ny, &y0out[i-1, 0], &k1[0], &k2[0], &k3[0], &k4[0],
                    1.0, h/6, h/3, h/3, h/6, &y0out[i, 0])


def rk4(CppReactionDiffusion rd, y0, tout):
    """
    simple explicit, fixed step size, Runge Kutta 4th order integrator.
    Use for debugging/testing.
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=2] y0out = np.empty((tout.size, y0.size), dtype=np.float64)
    cdef cnp.ndarray[cnp.float64_t, ndim=2] y1out = np.empty((tout.size, y0.size), dtype=np.float64)
    _rk4(rd.thisptr, np.asarray(y0), np.asarray(tout), y0out, y1out)
    return y0out.reshape((tout.size, rd.N, rd.n)), y1out.reshape((tout.size, rd.N, rd.n))
