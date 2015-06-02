# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from future.builtins import *

from itertools import product

import pytest

import simple

TR_FLS = (True, False)


@pytest.mark.parametrize('params', list(product(TR_FLS, TR_FLS)))
def test_simple(params):
    logy, logt = params
    simple.main(logy=logy, logt=logt)
