#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from four_species import main

if __name__ == '__main__':
    for implementation in ['py', 'cy', 'cpp']:
        print("Testing implementation={}".format(implementation))
        sys, y0, t0, ref_f, ref_J = main(implementation=implementation, demo=False, N=1)
        ny = sys.n*sys.N
        J = np.zeros((ny, ny))
        fout = np.zeros(ny)
        sys.f(0.0, y0, fout)
        assert np.allclose(fout, ref_f)
        sys.dense_jac_rmaj(0.0, y0, J, 1.0, False)

        np.set_printoptions(precision=3, linewidth=120)
        print "Calc J:"
        print J
        print "Ref J:"
        print ref_J

        assert np.allclose(J, ref_J)
