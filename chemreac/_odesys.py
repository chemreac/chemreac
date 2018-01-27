# -*- coding: utf-8 -*-

from pyodesys import ODESys as _ODESys
from pyodesys.results import Result
from .integrate import run


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
        pass
