# -*- coding: utf-8 -*-

from pyodesys import ODESys as _ODESys
from pyodesys.result import Result
from .integration import run

class ODESys(_ODESys):

    def __init__(self, rd):
        self.rd = rd

    names = property(lambda self: self.substance_names)
    latex_names = property(lambda self: self.substance_latex_names)

    def integrate(self, x, y0, params=None):
        integration = run(self.rd, y0, x, **kwargs)
        return Result(integr.tout, integr.Cout, params, integr.info, self)
