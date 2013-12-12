# -*- coding: utf-8 -*-

import numpy as np

def spy(sys):
    # Spy
    from matplotlib import pyplot as plt
    b = sys.spy()
    plt.spy(b)
    plt.show()

def coloured_log_spy(sys, t, y, h, cmap_name='gray'):
    import matplotlib.cm
    from matplotlib import pyplot as plt
    J = np.zeros([sys.n*sys.N]*2)
    sys.dense_jac_rmaj(t, y, J, factor=h, sub_one=True)
    plt.imshow(np.log(J), cmap=matplotlib.cm.get_cmap(cmap_name), interpolation='none')
    plt.colorbar()
    plt.show()

def wrap(cls, **wrap_kwargs):
    """
    Lazy overriding of keyword arguments to class initializing
    """
    def callback(*args, **kwargs):
        new_kwargs = kwargs.copy()
        new_kwargs.update(wrap_kwargs)
        return cls(*args, **new_kwargs)
    return callback
