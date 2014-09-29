
.PHONY: images-png

IMAGES=_generated/analytic_diffusion.png

_generated/analytic_diffusion.png: examples/examples/analytic_diffusion.py
	python $< --plot --nstencil 3 --k 0.1 --geom f --savefig $@

images-png: $(IMAGES)
