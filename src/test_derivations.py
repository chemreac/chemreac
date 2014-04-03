#!/usr/bin/env python
# -*- coding: utf-8 -*-

from sympy import *

from derivations import finite_diff

def test_finite_diff():
    x = IndexedBase('x')
    y = IndexedBase('y')
    i = Idx('i')

    # y = x**2
    dx_ = {i-1: 1, i: 1}.get
    dy_ = {i-1:-1, i: 1}.get

    xsubs = {x[i-1]: -1, x[i]: 0, x[i+1]: 1}

    Dy, DDy = finite_diff(x, y, i, 0, dx_, dy_)
    assert Dy == 0
    assert DDy == 2

    Dy, DDy = finite_diff(x, y, i, 1, dx_, dy_)
    assert (Dy - 2*(x[i+1]-x[i])).subs(xsubs) == 0
    assert DDy == 2

    Dy, DDy = finite_diff(x, y, i, -1, dx_, dy_)
    assert (Dy.subs(xsubs) + 2*(x[i]-x[i-1])).subs(xsubs) == 0
    assert DDy == 2

    # y = 3*x**2 + 2*x + 1, Dy = 6*x + 2, DDy = 6
    x_ = [2, 5, 9]
    y_ = lambda x: 3*x**2 + 2*x + 1
    Dy_ = lambda x: 6*x + 2
    DDy_ = lambda x: 6

    dx_ = {i-1: 3, i: 4}.get
    dy_ = {i-1: y_(x_[1])-y_(x_[0]), i: y_(x_[2])-y_(x_[1])}.get

    xsubs = dict(zip([x[i+j] for j in (-1,0,1)], x_))

    Dy, DDy = finite_diff(x, y, i, 0, dx_, dy_)
    assert Dy == Dy_(x_[1])
    assert DDy == DDy_(x_[1])

    Dy, DDy = finite_diff(x, y, i, 1, dx_, dy_)
    assert Dy.subs(xsubs) == Dy_(x_[2])
    assert DDy.subs(xsubs) == DDy_(x_[2])

    Dy, DDy = finite_diff(x, y, i, -1, dx_, dy_)
    assert Dy.subs(xsubs) == Dy_(x_[0])
    assert DDy.subs(xsubs) == DDy_(x_[0])


if __name__ == '__main__':
   test_finite_diff()
