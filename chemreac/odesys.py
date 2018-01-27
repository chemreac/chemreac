# -*- coding: utf-8 -*-

from pyodesys import ODESys as _ODESys
from pyodesys.results import Result
from .integrate import run


class ODESys(_ODESys):

    def __init__(self, rd):
        if rd.N > 1:
            raise NotImplementedError("ODESys expects single bin for now")
        self.rd = rd

    names = property(lambda self: self.rd.substance_names)
    latex_names = property(lambda self: self.rd.substance_latex_names)
    param_names = property(lambda self: self.rd.param_names)

    def integrate(self, x, y0, params=None, **kwargs):
        if isinstance(params, dict):
            params = [params[k] for k in self.param_names]
        if params is not None:
            self.rd.k = params
        integr = run(self.rd, [y0[k] for k in self.names], x, **kwargs)
        return Result(integr.tout, integr.Cout[:, 0, :], rd.k, integr.info, self)

    def chained_parameter_variation(self, durations, y0, varied_params, default_params=None,
                                    integrate_kwargs=None, x0=None, npoints=1, numpy=None):
        if list(varied_params) != ['doserate']:
            raise NotImplementedError("For now only varied doserate is supported")
        pass
