
.PHONY: images-png

IMAGES=_generated/analytic_diffusion.png _generated/steady_state.png _generated/steady_state.html

_generated/steady_state.png: examples/examples/steady_state.py
	python $< --plot --savefig $@

_generated/steady_state.html: examples/examples/steady_state.py
	python $< --plot --savefig $@

_generated/analytic_diffusion.png: examples/examples/analytic_diffusion.py
	python $< --plot --nstencil 3 --k 0.1 --geom f --savefig $@

images-png: $(IMAGES)
