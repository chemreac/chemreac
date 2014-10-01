
.PHONY: images-png

IMAGES=_generated/analytic_diffusion.png _generated/steady_state_approx.png _generated/steady_state_approx.html

_generated/steady_state_approx.png: examples/examples/steady_state_approx.py
	python $< --plot --savefig $@

_generated/steady_state_approx.html: examples/examples/steady_state_approx.py
	python $< --plot --savefig $@

_generated/analytic_diffusion.png: examples/examples/analytic_diffusion.py
	python $< --plot --nstencil 3 --k 0.1 --geom f --savefig $@

images-png: $(IMAGES)
