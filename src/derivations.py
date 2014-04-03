# -*- coding: utf-8 -*-

from sympy import *


def finite_diff(x, y, i, offset=0, dx_=None, dy_=None):
    """
    Returns finite difference formulae for first and second derivative
    of of y wrt x from three arbitrarily spaced points (parabola interpolation):
        [(x[i-1], y[i-1]), (x[i], y[i]), (x[i+1], y[i+1])]
    either at x[i+offset]:
      the left point (offset=-1) or
      the central point (offset=0) or
      the right point (offset=1)

    the separations between the points are: x[i]-x[i-1], x[i+1]-x[i]
    or if a `dx_` forward delta x callback (optional) is provided it
    will be called with `i-1` and `i` as arguments. In the same manner
    the user may pass `dy_` callback. This may be useful if furter symbolic
    tranformations are to be performed before
    """
    dx_ = dx_ or (lambda idx: x[idx+1]-x[idx])
    dy_ = dy_ or (lambda idx: y[idx+1]-y[idx])
    alpha, beta, a, b, c, d = symbols('alpha beta a b c d')
    # offset is -1 (left boundary), 0 anywhere inbetween, 1 (right boundary)
    coeff_mtx = Matrix([[beta**2, beta], [alpha**2, -alpha]])
    sol = coeff_mtx.inv().dot([c,d])
    subs = {alpha: dx_(i-1), beta: dx_(i), c: dy_(i), d: -dy_(i-1)}
    a_ = sol[0].subs(subs).simplify()
    b_ = sol[1].subs(subs).simplify()
    return 2*a_*{-1: -dx_(i-1), 0: 0, 1: dx_(i)}.get(offset) + b_, 2*a_


t, r, D = symbols('t r D', real=True, positive=True)
phi = Function('phi')(t,r)
phi_t = phi.diff(t)
phi_r = phi.diff(r,1)
phi_rr = phi.diff(r,2)
diffusion_eq = {
    'FLAT': Eq(phi_t, D*phi_rr),
    'CYLINDRICAL': Eq(phi_t, D*(phi_rr+1/r*phi_r)),
    'SPHERICAL': Eq(phi_t, D*(phi_rr+2/r*phi_r)),
}
geoms = ['FLAT', 'CYLINDRICAL', 'SPHERICAL']

N = Symbol('N', integer=True)
i = Idx('i', N) # 0 <= i <= N-1
x = IndexedBase('xc', shape=(N,)) # x is boundaries in code, xc is bin-centers
dx, dx_ = IndexedBase('delta', shape=(N,)), lambda idx: x[idx+1]-x[idx]
y = IndexedBase('y')
gamma = IndexedBase(Symbol('gamma'))
gamma_ = lambda idx: y[idx+1]-y[idx]
delta_ = lambda idx: dx[idx]


def write_code(offset, kind, dct):
    Dy_, DDy_ = finite_diff(x, y, i, offset, delta_, gamma_)
    # xsymb = Symbol('x')
    # yfunc = Function('y')(xsymb)
    # Dy, DDy = yfunc.diff(xsymb), yfunc.diff(xsymb,2)
    Delta, Delta_ = IndexedBase('Delta'), lambda idx: dx[i-1]+dx[i] #x[i+1]-x[i-1]
    yd = y_m, y_n, y_p = symbols('y_m y_n y_p')
    ys, dummies = (y[i-1], y[i], y[i+1]), yd
    fw_subs = dict(zip(ys, dummies))
    bw_subs = dict(zip(dummies, ys))
    canonical = {dx[i]: dx_(i), dx[i-1]: dx_(i-1)}
    canonical[Delta[i]] = Delta_(i).subs(canonical)

    expr = diffusion_eq[kind].subs({phi_r: Dy_, phi_rr: DDy_, r: x[i]})
    expr = expr.subs({Delta_(i): Delta[i]})
    expr = Eq(expr.lhs/D, expr.rhs/D)
    fexpr = expr.subs({gamma_(i): gamma[i], gamma_(i-1): gamma[i-1]})
    fexpr = fexpr.collect([gamma[i], gamma[i-1]])
    #display(fexpr)#,gamma[i-1]]))
    f_in_gamma = {
        gamma[i]: fexpr.rhs.expand().coeff(gamma[i]),
        gamma[i-1]: fexpr.rhs.expand().coeff(gamma[i-1])
    }
    for j, gamma_idx in enumerate([i-1,i]):
        code = ccode(f_in_gamma[gamma[gamma_idx]].subs(
            canonical), contract=False)
        dct[kind, offset, ('bw', 'fw')[j]] = code
    from operator import add as add_
    assert (fexpr.rhs - reduce(add_, [k*v for k,v in f_in_gamma.items()])).simplify() == 0
    jac=[]
    for j, dy in enumerate(yd):
        Jexpr=expr.rhs.subs(fw_subs).diff(dy).simplify().subs(bw_subs)#.subs(canonical)
        canonical_Jexpr = Jexpr.subs(canonical)
        code = ccode(canonical_Jexpr, contract=False)
        dct[kind, offset, ('Jp','Jc','Jn')[j]] = code
        jac.append(code)

def main():
    import itertools
    d={}
    for offset, kind in itertools.product([-1,0,1], geoms):
        print(offset, kind)
        write_code(offset, kind, d)
    import pickle
    pickle.dump({'COEFF': d}, open('exprs.pkl', 'wb'))
    cmdstr = '''python -c "import pickle; print(pickle.load(open('exprs.pkl', 'rb'))['flat',-1,'bw'])"'''
    print("Run e.g.: "+cmdstr)

if __name__ == '__main__':
   main()
