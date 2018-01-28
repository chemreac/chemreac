# -*- coding: utf-8 -*-
import time
import numpy as np
from pyodesys import ODESys as _ODESys
from pyodesys.results import Result
from .integrate import run
from ._chemreac import cvode_predefined_durations_modulations


class ODESys(_ODESys):

    def __init__(self, rd, k_from_params=None):
        if rd.N > 1:
            raise NotImplementedError("ODESys expects single bin for now")
        self.rd = rd
        self.k_from_params = k_from_params

    names = property(lambda self: self.rd.substance_names)
    latex_names = property(lambda self: self.rd.substance_latex_names)
    param_names = property(lambda self: self.rd.param_names)
    # dep_by_name = True
    # par_by_name = True

    def integrate(self, x, y0, params=None, **kwargs):
        if params is not None:
            self.rd.k = self.k_from_params(params)
        integr = run(self.rd, [y0[k] for k in self.names], x, **kwargs)
        return Result(integr.tout, integr.Cout[:, 0, :], None, integr.info, self)

    def chained_parameter_variation(self, durations, y0, varied_params, default_params=None,
                                    integrate_kwargs=None, x0=None, npoints=1, numpy=None):
        if list(varied_params) != ['doserate']:
            raise NotImplementedError("For now only varied doserate is supported")
        if default_params is not None:
            self.rd.k = self.k_from_params(default_params)
        if x0 is not None:
            assert x0 == 0
        integrate_kwargs = integrate_kwargs or {}
        atol = integrate_kwargs.pop('atol', 1e-8)
        if isinstance(atol, float):
            atol = [atol]
        rtol = integrate_kwargs.pop('rtol', 1e-8)
        method = integrate_kwargs.pop('method', 'bdf')
        time_cpu = time.clock()
        time_wall = time.time()
        tout, yout = cvode_predefined_durations_modulations(
            self.rd, np.array([y0[k] for k in self.names], dtype=np.float64),
            np.array(durations, dtype=np.float64),
            np.array(varied_params['doserate'], dtype=np.float64),
            atol=atol, rtol=rtol, method=method, npoints=npoints, **integrate_kwargs)
        info = dict(
            nfev=self.rd.nfev,
            njev=self.rd.njev,
            time_wall=time.time() - time_wall,
            time_cpu=time.clock() - time_cpu,
            success=True
        )
        return Result(tout, yout[:, 0, :], None, info, self)
