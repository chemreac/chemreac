# -*- coding: utf-8 -*-

import numpy as np

def spy(sys):
    # Spy
    from matplotlib import pyplot as plt
    b = sys.spy()
    plt.spy(b)
    plt.show()

def coloured_spy(A, cmap_name='gray', ax = None):
    from matplotlib import pyplot as plt
    from matplotlib.ticker import MaxNLocator
    from matplotlib.cm import get_cmap

    if not ax:
        # fig = plt.figure()
        # ax = fig.add_subplot(111)
        plt.imshow(A, cmap=get_cmap(cmap_name), interpolation='none')
        ax = plt.gca()

    ya = ax.get_yaxis()
    ya.set_major_locator(MaxNLocator(integer=True))
    xa = ax.get_xaxis()
    xa.set_major_locator(MaxNLocator(integer=True))
    plt.colorbar()
    return ax
