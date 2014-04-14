#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division, absolute_import

import argh
import numpy as np
import matplotlib.pyplot as plt

from chemreac import ReactionDiffusion, FLAT, CYLINDRICAL, SPHERICAL, Geom_names
from chemreac.integrate import run

def flat_analytic(x, t, D, mu):
    return (4*np.pi*D*t)**-0.5 * np.exp(-(x-mu)**2/(4*D*t))

def spherical_analytic(x, t, D, mu):
    return (4*np.pi*D)**-0.5 * t**-1.5 * np.exp(-(x-mu)**2/(4*D*t))

def cylindrical_analytic(x, t, D, mu):
    return (4*np.pi*D*t)**-1 * np.exp(-(x-mu)**2/(4*D*t))


def interleave(arrays, axis=0):
    """
    Interleaves multiple arrays along a user specified axis.

    Parrameters
    -----------
    arrays : sequence of array_like with same shape and dtype
         Input arrays, if a sequence of length 1 is given, the fist
    axis : int
         Along what axis interleaving should be performed

    Returns
    -------
    interleave: ndarray
        Interleaved array of rank queal to arrays in `arrays` with shape
        in `axis` dimension multiplied with number of arrays in `arrays`.

    Raises
    ------
    ValueError
        if axis is larger than ndim of arrays

    See Also
    --------
    concatenate

    Examples
    --------
    >>> a = np.arange(3)
    >>> interleave(a, a+6)
    array([0, 6, 1, 7, 2, 8])

    >>> a = np.arange(4).reshape((2,2))
    >>> interleave(a, a+4, axis=1)
    array([[0, 4, 1, 5],
           [2, 6, 3, 7]])
    """
    narrays = len(arrays)
    if narrays == 0:
        raise ValueError("interleave takes at least one (reasonably two) arrays")
    arrays = map(np.asarray, arrays)
    if any(x.dtype != arrays[0].dtype for x in arrays[1:]):
        raise ValueError("dtype of arrays must be consistent")
    if any(x.shape != arrays[0].shape for x in arrays[1:]):
        raise ValueError("Shapes of arrays must be consistent")
    nd = arrays[0].ndim
    if nd == 0: raise ValueError("interleave only works on arrays of 1 or more dimensions")
    if axis < 0:
        axis += nd
    if (axis >= nd):
        raise ValueError("axis must be less than arr.ndim; axis=%d, rank=%d."
            % (axis, nd))

    c = np.empty([narrays*n if i==axis else n for i, n \
                  in enumerate(arrays[0].shape)], dtype=arrays[0].dtype)
    for i, arr in enumerate(arrays):
        slicing = [slice(i, arrays[0].shape[j]*narrays, narrays) if\
                   j==axis else Ellipsis for j in range(arrays[0].ndim)]
        c[slicing] = arr
    return c


def main(D=2e-3, t0=3., tend=7., x0=0., xend=1., mu=None, N=2048, nt=30, geom='f',
         logt=False, logy=False, random=False, k=0.0):
    decay = (k != 0.0)
    mu = float(mu or x0)
    tout = np.linspace(t0, tend, nt)

    assert geom in 'fcs'
    geom = {'f': FLAT, 'c': CYLINDRICAL, 's': SPHERICAL}[geom]
    print(Geom_names[geom])
    analytic = {
        FLAT: flat_analytic,
        CYLINDRICAL: cylindrical_analytic,
        SPHERICAL: spherical_analytic
    }[geom]

    # Steup the system
    x = np.linspace(x0, xend, N+1)
    if random: x += (np.random.random(N+1)-0.5)*(xend-x0)/(N+2)
    sys = ReactionDiffusion(
        2 if decay else 1,
        [[0]] if decay else [],
        [[1]] if decay else [],
        [k] if decay else [],
        N,
        D=[D]*2 if decay else [D],
        x=x, geom=geom, logy=logy, logt=logt)

    # Calc initial conditions / analytic reference values
    t = tout.copy().reshape((nt,1))
    yref = analytic(sys.xc, t, D, mu)
    if decay: yref = interleave((yref*np.exp(-k*t), yref*(1-np.exp(-k*t))), axis=1)
    y0 = yref[0,:]

    # Run the integration
    y = np.log(y0) if logy else y0
    t = np.log(tout) if logt else tout
    yout, info = run(sys, y, t, atol=1e-6, rtol=1e-6, with_jacobian=False, method='bdf')
    if logy: yout = np.exp(yout)
    print(info)

    # Plot results
    def plot(y, c, ttl=None):
        plt.plot(sys.xc, y, c=c)
        if N < 100: plt.vlines(sys.x, 0, np.ones_like(sys.x)*y0[0], linewidth=.1, colors='gray')
        plt.xlabel('x / m')
        plt.ylabel('C / M')
        if ttl: plt.title(ttl)

    stride = 2 if decay else 1
    for i in range(nt):
        c = 1-tout[i]/tend
        c = (1.0-c, .5-c/2, .5-c/2)

        plt.subplot(4,1,1)
        plot(yout[i,::stride], c, 'Simulation (N={})'.format(sys.N))
        if decay: plot(yout[i,1::stride], c[::-1])

        plt.subplot(4,1,2)
        plot(yref[i,::stride], c, 'Analytic')
        if decay: plot(yref[i,1::stride], c[::-1])

        plt.subplot(4,1,3)
        plot((yref[i,::stride]-yout[i,::stride])/info['atol'], c,
             'Abs. err. / Abs. tol. (={})'.format(info['atol']))
        if decay: plot((yref[i,1::stride]-yout[i,1::stride])/info['atol'], c[::-1])

    plt.subplot(4,1,4)
    plt.plot(tout, np.sum((yref[:,::stride]-yout[:,::stride])**2/N, axis=1)**0.5/info['atol'])
    if decay: plt.plot(tout, np.sum((yref[:,::stride]-yout[:,1::stride])**2/N,
                                    axis=1)**0.5/info['atol'])
    plt.xlabel('Time / s')
    plt.ylabel(r'$\sqrt{\langle E^2 \rangle} / atol$')
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    argh.dispatch_command(main)
