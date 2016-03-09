
.PHONY: images-png

IMAGES=_generated/analytic_diffusion.png _generated/steady_state_approx.png _generated/steady_state_approx.html _generated/aqueous_radiolysis.png _generated/analytic_N_scaling.png _generated/auto_efield.png _generated/decay.png _generated/decay_long.png _generated/decay_long_damp.png _generated/four_species.png _generated/equilibrium.png _generated/equilibrium_unscaled.png _generated/equilibrium_scaled.png _generated/const_surf_conc.png _generated/const_surf_conc_logy_logx.png

_generated/steady_state_approx.png: examples/examples/steady_state_approx.py
	python $< --plot --savefig $@

_generated/steady_state_approx.html: examples/examples/steady_state_approx.py
	python $< --plot --savefig $@

_generated/analytic_diffusion.png: examples/examples/analytic_diffusion.py
	python $< --plot --nstencil 3 --nspecies 3 --geom f --savefig $@

_generated/aqueous_radiolysis.png: examples/examples/aqueous_radiolysis.py
	python $< --doserate 25 --plot --savefig $@

_generated/analytic_N_scaling.png: examples/examples/analytic_N_scaling.py
	python $< --nNs 6 --plot --savefig $@

_generated/auto_efield.png: examples/examples/auto_efield.py
	python $< --plot --savefig $@

_generated/decay.png: examples/examples/decay.py
	python $< --plot --savefig $@

_generated/decay_long.png: examples/examples/decay.py
	python $< --plot --savefig $@ --rates 1.0 --logy --logt --rtol 1e-13 --atol 1e-6 --scale-err 100.0 --plotlogy --nt 1024 --tend 1700

_generated/decay_long_damp.png: examples/examples/decay.py
	python $< --plot --savefig $@ --rates 1.0 --logy --logt --rtol 1e-13 --atol 1e-6 --scale-err 100.0 --plotlogy --nt 1024 --tend 1700 --sigm-damp

_generated/four_species.png: examples/examples/four_species.py
	python $< --plot --savefig $@

_generated/equilibrium.png: examples/examples/equilibrium.py
	python $< --plot --savefig $@

_generated/equilibrium_unscaled.png: examples/examples/equilibrium.py
	python $< --A0 1.0 --B0 1e-10 --C0 1e-30 --kf 10 --kb 1 --t0 0 --tend 5 --plot --plotlogy --plotlogt --savefig $@

_generated/equilibrium_scaled.png: examples/examples/equilibrium.py
	python $<  --scaling 1e10 --A0 1.0 --B0 1e-10 --C0 1e-30 --kf 10 --kb 1 --t0 0 --tend 5 --plot --plotlogy --plotlogt --savefig $@

_generated/const_surf_conc.png: examples/examples/const_surf_conc.py
	python $< --plot --savefig $@

_generated/const_surf_conc_logy_logx.png: examples/examples/const_surf_conc.py
	python $< --logx --logy --x0 1e-6 --scaling 1e-20 --factor 1e12 --plot --savefig $@


images-png: $(IMAGES)
