// This is a templated source file. (Make sure you're editing the 
// template)
// Render template using Mako (Python templating engine)

#include <algorithm> // count
#include <vector>
#include "cpp_chem.hpp"

#ifdef DEBUG
#include <cstdio>
#endif

%if USE_OPENMP:
#include <omp.h>
#ifndef _OPENMP
  #error
#endif
%else:
#define omp_get_thread_num() 0
#ifdef _OPENMP
  #error
#endif
%endif

using std::vector;
using std::count;

// 1D discretized reaction diffusion
ReactionDiffusion::ReactionDiffusion(
    int n,
    int N, 
    vector<vector<int> > stoich_reac,
    vector<vector<int> > stoich_prod,
    vector<vector<int> > stoich_actv,
    vector<double> k,
    vector<double> D,
    vector<double> x, // separation
    int mode,
    int geom_
    ):
    n(n), N(N), stoich_reac(stoich_reac), 
    stoich_prod(stoich_prod), stoich_actv(stoich_actv),
    k(k), D(D), x(x), mode(mode)
{
    if (stoich_reac.size() != stoich_prod.size())
        throw std::length_error(
            "stoich_reac and stoich_prod of different sizes.");
    if (k.size() != stoich_prod.size())
        throw std::length_error(
            "k and stoich_prod of different sizes.");
    if (N>1){
        if (D.size() != n)
            throw std::length_error(
                "Length of D does not match number of species.");
        if (x.size() != N + 1)
            throw std::length_error(
                "Number bin edges != number of compartments + 1.");
    }

    switch(geom_) {
    case 0: geom = Geom::FLAT; break;
    case 1: geom = Geom::SPHERICAL; break;
    case 2: geom = Geom::CYLINDRICAL; break;
    }

    dx = new double[N-1];

    for (int i=0; i<N-1; ++i)
        dx[i] = (x[i+2]-x[i])/2;

    nr = stoich_reac.size();

    coeff_reac = new int[nr*n];
    coeff_prod = new int[nr*n];
    coeff_totl = new int[nr*n];
  
    for (int ri=0; ri<nr; ++ri){ // reaction index 
        if (stoich_actv[ri].size() == 0)
            stoich_actv[ri] = stoich_reac[ri]; // Strict massaction
        for (int si=0; si<n; ++si){ // species index
            coeff_reac[ri*n+si] = count(stoich_reac[ri].begin(), 
                                        stoich_reac[ri].end(), si);
            coeff_prod[ri*n+si] = count(stoich_prod[ri].begin(), 
                                        stoich_prod[ri].end(), si);
            coeff_totl[ri*n+si] = coeff_prod[ri*n+si] -\
                coeff_reac[ri*n+si];
        }
    }
    nfeval = 0;
    njeval = 0;
}

ReactionDiffusion::~ReactionDiffusion()
{
    delete []dx;
    delete []coeff_reac;
    delete []coeff_prod;
    delete []coeff_totl;
}

void
ReactionDiffusion::_fill_local_r(const double * const restrict yi, double * const restrict local_r)
{
    // intent(out) :: local_r
    for (int j=0; j<nr; ++j){
        // reaction j
        local_r[j] = k[j];
        for (int l=0; l<stoich_actv[j].size(); ++l){
            // reactant index l
            int m = stoich_actv[j][l];
            local_r[j] *= yi[m];
        }
    }
}


double
ReactionDiffusion::flux(int bi, int si, const double * const restrict y)
{
    // bi: bin index, si: species index
    return (y[(bi+1)*n+si] - y[bi*n+si])/dx[bi];
}

// 4π cancel
#define A(i) (x[i]*x[i]) // 4πrₖ²
#define V(i) ((x[i+1]*x[i+1]*x[i+1] - x[i]*x[i]*x[i])/3) // 4πrₖ₊₁³/3 - 4πrₖ³/3
double
ReactionDiffusion::diffusion_contrib(int bi, int si, const double * const restrict fluxes)
{
    // bi: bin index, si: species index, fluxes: mol/m2/s through right wall of bin
    const double * const flux = fluxes + bi*n + si;
    double contrib = 0;
    switch(geom){
    case Geom::FLAT :
        if (bi > 0)
            contrib += fluxes[(bi-1)*n+si]/dx[bi];
        if (bi < N-1)
            contrib -= fluxes[bi*n+si]/dx[bi];
        break;
    case Geom::SPHERICAL :
        if (bi > 0)
            contrib += flux[(bi-1)*n+si]*A(bi)/V(bi);
        if (bi < N-1)
            contrib -= flux[bi*n+si]*A(bi+1)/V(bi);
        break;
    }
    
    return D[si]*contrib;
}

double
ReactionDiffusion::diffusion_contrib_jac_prev(
    int bi, int si, const double * const restrict y)
{
    switch(geom){
    case Geom::FLAT :
        return -D[si]/(dx[bi-1]*dx[bi]);
    case Geom::SPHERICAL :
        return -D[si]/dx[bi-1]*A(bi)/V(bi);
    }
}

double
ReactionDiffusion::diffusion_contrib_jac_next(
    int bi, int si, const double * const restrict y)
{
    switch(geom){
    case Geom::FLAT :
        return -D[si]/(dx[bi]*dx[bi]);
    case Geom::SPHERICAL :
        return -D[si]/dx[bi]*A(bi+1)/V(bi);
    }
}
#undef A
#undef V


void
ReactionDiffusion::f(double t, const double * const restrict y, double * const restrict dydt)
{
    // initialize fluxes
    double * const fluxes = new double[(N-1)*n];
    ${"#pragma omp parallel for" if USE_OPENMP else ""}
    for (int i=0 i<N-1; ++i)
        for (int k=0; k<n; ++k)
            fluxes[i*n + k] = flux(i, k, y);

    ${"double * local_r = new double[nr];" if not USE_OPENMP else ""}
    ${"#pragma omp parallel for" if USE_OPENMP else ""}
    for (int i=0; i<N; ++i){
        // compartment i
        ${"double * local_r = new double[nr];" if USE_OPENMP else ""}

        for (int j=0; j<n; ++j)
            dydt[i*n+j] = 0.0; // zero out

        // Contributions from reactions
        // ----------------------------
        _fill_local_r(y+i*n, local_r);
        for (int j=0; j<nr; ++j){
            // reaction j
            for (int k=0; k<n; ++k){
                // species k
                int overall = coeff_totl[j*n + k];
                if (overall != 0)
                    dydt[i*n + k] += overall*local_r[j];
            }
        }

        if (N>1){
            // Contributions from diffusion
            // ----------------------------
            for (int si=0; si<n; ++si)// species index si
                dydt[i*n + si] += diffusion_contrib(i, si, fluxes);
        }

        ${"delete []local_r;" if USE_OPENMP else ""}

    }

    ${"delete []local_r;" if not USE_OPENMP else ""}
    delete []fluxes;
    nfeval++;
}


%for token, imaj, imin in                                       \
    [('dense_jac_rmaj', '((bri)*n + ri)', '(bci)*n + ci'), \
     ('dense_jac_cmaj', '((bci)*n + ci)', '(bri)*n + ri'),\
     ('banded_jac_cmaj','((bci)*n + ci)', '2*n+((bri)*n+ri)-((bci)*n+ci)'),\
     ('banded_packed_jac_cmaj','((bci)*n + ci)', 'n+((bri)*n+ri)-((bci)*n+ci)'),\
        ]:
#define JAC(bri, bci, ri, ci) ja[${imaj}*ldj+${imin}]
void
ReactionDiffusion::${token}(double t, const double * const restrict y,
                            double * const restrict ja, int ldj)
{
    // intent(out) :: ja

    // initialize fluxes
    double * const fluxes = new double[(N-1)*n];
    ${"#pragma omp parallel for" if USE_OPENMP else ""}
    for (int i=0 i<N-1; ++i)
        for (int k=0; k<n; ++k)
            fluxes[i*n + k] = flux(i, k, y);

#ifndef _OPENMP
    double * local_r = new double[nr];
#else
#pragma omp parallel for
#endif
    for (int bi=0; bi<N; ++bi){
        const double * const C = y + bi*n; // Conc. in `bi:th` compartment

#ifdef _OPENMP
        double * local_r = new double[nr];
#endif
    
        // Contributions from reactions
        // ----------------------------
        _fill_local_r(C, local_r);
        for (int j=0; j<n; ++j){
            // species j
            for (int m=0; m<n; ++m){
                // derivative wrt species m
                // j_i[j, m] = Sum_l(n_lj*Derivative(r[l], C[m]))
                JAC(bi,bi,j,m) = 0.0;
                for (int l=0; l<nr; ++l){
                    // reaction l
                    if (coeff_totl[l*n + j] == 0)
                        continue;
                    else if (coeff_reac[l*n + m] == 0)
                        continue;
                    JAC(bi,bi,j,m) += coeff_totl[l*n + j]*\
                        coeff_reac[l*n + m]*local_r[l]/C[m];
                }
            }
        }

        if (N>1){
            // Contributions from diffusion
            // ----------------------------
            for (int si=0; si<n; ++si){
                // species index si
                if (bi > 0){
                    double tmp = diffusion_contrib_jac_prev(bi, si, y);
                    JAC(bi, bi-1, si, si)  = tmp;
                    JAC(bi, bi,   si, si) -= tmp; // from symmetry
                }
                if (bi < N-1){
                    double tmp = diffusion_contrib_jac_next(bi, si, y);
                    JAC(bi, bi+1, si, si)  = tmp;
                    JAC(bi, bi,   si, si) += tmp; // from symmetry
                }
            }
        }
#ifdef _OPENMP
        delete []local_r;
#endif
    }

#ifndef _OPENMP
    delete []local_r;
#endif

    delete []fluxes;
    njeval++;
}
#undef JAC
%endfor
