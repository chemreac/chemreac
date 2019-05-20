# -*- coding: utf-8 -*-
import time
import numpy as np
from pyodesys import ODESys as _ODESys
from pyodesys.results import Result
from chempy.units import get_derived_unit, unitless_in_registry, uniform, patched_numpy as pnp
from .integrate import run
from ._chemreac import cvode_predefined_durations_fields


class ODESys(_ODESys):

    def __init__(self, rd, k_from_params=None, variables_from_params=None):
        if rd.N > 1:
            raise NotImplementedError("ODESys expects single bin for now")
        self.rd = rd
        self.k_from_params = k_from_params
        self.variables_from_params = variables_from_params

    ny = property(lambda self: self.rd.n*self.rd.N)
    names = property(lambda self: self.rd.substance_names)
    latex_names = property(lambda self: self.rd.substance_latex_names)
    param_names = property(lambda self: self.rd.param_names)
    autonomous_interface = property(lambda self: not self.rd.logt)
    numpy = pnp
    # dep_by_name = True
    # par_by_name = True

    def _get_units_util(self):
        if self.rd.unit_registry is None:
            _dedim = lambda x: np.array(x)
            time_u = 1
            conc_u = 1
            dr_u = 1
        else:
            _dedim = lambda x: unitless_in_registry(x, self.rd.unit_registry)
            time_u = get_derived_unit(self.rd.unit_registry, 'time')
            conc_u = get_derived_unit(self.rd.unit_registry, 'concentration')
            dr_u = get_derived_unit(self.rd.unit_registry, 'doserate')
        return locals()

    def integrate(self, x, y0, params=None, integrator='cvode', **kwargs):
        if params is not None and self.k_from_params is not None:
            self.rd.k = self.k_from_params(self, params)
        if 'doserate' in (params or {}):
            self.rd.set_with_units(
                'fields', [[self.variables_from_params['density'](self, params)*params['doserate']]])
        if 'atol' in kwargs and isinstance(kwargs['atol'], dict):
            kwargs['atol'] = [kwargs['atol'][k] for k in self.names]
        integr = run(self.rd, [y0[k] for k in self.names] if isinstance(y0, dict) else y0,
                     x, integrator=integrator, **kwargs)
        pout = [params[k] for k in self.param_names] if self.param_names else None
        return Result(integr.with_units('tout'), integr.with_units('Cout')[:, 0, :],
                      pout, integr.info, self)

    def chained_parameter_variation(self, durations, y0, varied_params, default_params=None,
                                    integrate_kwargs=None, x0=None, npoints=1, numpy=None):
        if list(varied_params) != ['doserate']:
            raise NotImplementedError("For now only varied doserate is supported")
        if self.param_names != ['doserate']:
            raise NotImplementedError("We expect doserate to be varied for now")
        uutil = self._get_units_util()
        _dedim, time_u, conc_u, dr_u = [uutil[k] for k in '_dedim time_u conc_u dr_u'.split()]
        density = _dedim(self.variables_from_params['density'](self, default_params))
        if default_params:
            self.rd.k = _dedim(self.k_from_params(self, default_params))

        if x0 is not None:
            assert x0 == 0*time_u
        integrate_kwargs = integrate_kwargs or {}
        atol = integrate_kwargs.pop('atol', 1e-8)
        if isinstance(atol, float):
            atol = [atol]
        elif isinstance(atol, dict):
            atol = [atol[k] for k in self.names]
        rtol = integrate_kwargs.pop('rtol', 1e-8)
        method = integrate_kwargs.pop('method', 'bdf')
        integrator = integrate_kwargs.pop('integrator', 'cvode')
        if integrator != 'cvode':
            raise NotImplementedError("chained_parameter_variation requires cvode for now")
        drate = uniform(varied_params['doserate'])
        time_cpu = time.process_time()
        time_wall = time.time()
        tout, yout = cvode_predefined_durations_fields(
            self.rd, _dedim([y0[k] for k in self.names]),
            _dedim(durations),
            _dedim(drate*density),
            atol=atol, rtol=rtol, method=method, npoints=npoints, **integrate_kwargs)
        info = dict(
            nsteps=-1,
            nfev=self.rd.nfev,
            njev=self.rd.njev,
            time_wall=time.time() - time_wall,
            time_cpu=time.process_time() - time_cpu,
            success=True,
            integrator=[integrator],
            t0_set=False,
            linear_solver=0,  # pyodesys.results.Result work-around for now (not important)
        )
        info.update(self.rd.last_integration_info)
        dr_out = np.concatenate((np.repeat(drate, npoints), drate[-1:]))
        return Result(tout*time_u, yout[:, 0, :]*conc_u, dr_out.reshape((-1, 1))*dr_u, info, self)
