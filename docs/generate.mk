
.PHONY: images-png

IMAGES=_generated/analytic_diffusion.png _generated/steady_state_approx.png _generated/steady_state_approx.html _generated/aqueous_radiolysis.png _generated/analytic_N_scaling.png _generated/auto_efield.png _generated/decay.png _generated/four_species.png _generated/equilibrium.png _generated/equilibrium_unscaled.png _generated/equilibrium_scaled.png

_generated/steady_state_approx.png: examples/examples/steady_state_approx.py
	python $< --plot --savefig $@

_generated/steady_state_approx.html: examples/examples/steady_state_approx.py
	python $< --plot --savefig $@

_generated/analytic_diffusion.png: examples/examples/analytic_diffusion.py
	python $< --plot --nstencil 3 --k 0.1 --geom f --savefig $@

_generated/aqueous_radiolysis.png: examples/examples/aqueous_radiolysis.py
	python $< --doserate 25 --plot --savefig $@

_generated/analytic_N_scaling.png: examples/examples/analytic_N_scaling.py
	python $< --nNs 6 --plot --savefig $@

_generated/auto_efield.png: examples/examples/auto_efield.py
	python $< --plot --savefig $@

_generated/decay.png: examples/examples/decay.py
	python $< --plot --savefig $@

_generated/four_species.png: examples/examples/four_species.py
	python $< --plot --savefig $@

_generated/equilibrium.png: examples/examples/equilibrium.py
	python $< --plot --savefig $@

_generated/equilibrium_unscaled.png: examples/examples/equilibrium.py
	python $< --A0 1.0 --B0 1e-10 --C0 1e-30 --kf 10 --kb 1 --t0 0 --tend 5 --plot --plotlogy --plotlogt --savefig $@

_generated/equilibrium_scaled.png: examples/examples/equilibrium.py
	python $<  --scaling 1e10 --A0 1.0 --B0 1e-10 --C0 1e-30 --kf 10 --kb 1 --t0 0 --tend 5 --plot --plotlogy --plotlogt --savefig $@


images-png: $(IMAGES)
