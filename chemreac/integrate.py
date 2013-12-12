# -*- coding: utf-8 -*-

from scipy.integrate import ode


class ScipyIntegrator:

    order = 5 # It actually uses dynamic order

    atol = None
    rtol = None

    pickle_attrs = ('ny', 'tout', 'eout', 'yout')

    def __getstate__(self):
        return [getattr(self, x) for x in self.pickle_attrs]

    def __setstate__(self, state):
        for attr, val in zip(self.pickle_attrs, state):
            setattr(self, attr, val)

    def f_vec(self, t, y):
        out = np.empty(self.ny)
        self.sys.f(t, y, out)
        return out



    def __init__(self, sys, y0, t0, integrator='vode',
                 method='adams', **kwargs):
        self.sys = sys
        self.ny = len(y0)
        self.__dict__.update(kwargs)

        self._ode = ode(self.f_vec, jac=self.dense_jac_rmaj_scipy)
        # Set stricter tolerances (used as a reference)
        self._ode.set_integrator(
            integrator, method=method,
            with_jacobian=True,
            atol=1e-12, rtol=1e-12)
        self._ode.set_initial_value(np.array(y0), t0)

    def dense_jac_rmaj_scipy(self, t, y):
        """
        Returns a dense jacobian matrix
        """
        out = np.zeros((self.ny, self.ny))
        self.sys.dense_jac_rmaj(t, y, out, 1.0, False)
        return out


    def _run(self, tend, nt, h=None, hmax=None, errcalc_freq=10):
        if nt == 0:
            # Dynamic step-size
            self._ode._integrator.iwork[2] = -1
            tout, yout, fout, eout, hout = [], [], [], [], []
            warnings.filterwarnings("ignore", category=UserWarning)
            zero_arr = np.zeros(self.ny)
            while self._ode.t < tend:
                fvec = self.f_vec(self._ode.t, self._ode.y)
                yout.append([self._ode.y, fvec])
                tout.append(self._ode.t)
                # Step size last used
                hout.append(self._ode._integrator.rwork[11-1])
                # Doesn't work (not ciritical)
                eout.append(self._ode._integrator.rwork[-self.ny:])
                self._ode.integrate(tend, step=True)
            warnings.resetwarnings()

            self.tout = np.array(tout)
            self.yout = np.array(yout).swapaxes(0,1)
            self.eout = np.array(eout)
            self.hout = np.array(hout)
        else:
            fvec = self.f_vec(self.tout[0], self._ode.y)
            self.tout[0] = self.tout[0]
            self.yout[0, 0, :] = self._ode.y
            self.yout[1, 0, :] = fvec
            self.eout[0, :] = 0.0
            for i, t in np.ndenumerate(self.tout[1:]):
                self._ode.integrate(t)
                fvec = self.f_vec(self._ode.t, self._ode.y)
                self.tout[i[0]+1] = self._ode.t
                self.yout[0, i[0]+1, :] = self._ode.y
                self.yout[1, i[0]+1, :] = fvec
                if self.use_extrapol_jac:
                    self.yout[2, i[0]+1, :] = self.f2_vec(
                        self._ode.t, self._ode.y, fvec)
                else:
                    self.yout[2, i[0], :] = 0.0
                self.eout[i[0]+1, :] = self.atol
