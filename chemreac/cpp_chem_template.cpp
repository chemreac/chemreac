## -*- coding: utf-8 -*-
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
    vector<vector<int> > stoich_actv_,
    vector<double> k,
    vector<double> D,
    vector<double> x, // separation
    int mode,
    int geom_
    ):
    n(n), N(N), stoich_reac(stoich_reac), 
    stoich_prod(stoich_prod), k(k), D(D), x(x),
    mode(mode)
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
    default: throw std::logic_error("Unknown geom.");
    }

    dx = new double[N-1];

    for (int i=0; i<N-1; ++i)
        dx[i] = (x[i+2]-x[i])/2;

    nr = stoich_reac.size();

    coeff_reac = new int[nr*n];
    coeff_prod = new int[nr*n];
    coeff_totl = new int[nr*n];
    coeff_actv = new int[nr*n];


    stoich_actv.reserve(nr);
    for (int rxni=0; rxni<nr; ++rxni){ // reaction index 
        if (stoich_actv_[rxni].size() == 0)
            stoich_actv.push_back(stoich_reac[rxni]); // Strict massaction
	else
	    stoich_actv.push_back(stoich_actv_[rxni]);
        for (int si=0; si<n; ++si){ // species index
            coeff_reac[rxni*n+si] = count(stoich_reac[rxni].begin(), 
                                        stoich_reac[rxni].end(), si);
            coeff_actv[rxni*n+si] = count(stoich_actv[rxni].begin(), 
                                        stoich_actv[rxni].end(), si);
            coeff_prod[rxni*n+si] = count(stoich_prod[rxni].begin(), 
                                        stoich_prod[rxni].end(), si);
            coeff_totl[rxni*n+si] = coeff_prod[rxni*n+si] -\
                coeff_reac[rxni*n+si];
        }
    }
}

ReactionDiffusion::~ReactionDiffusion()
{
    delete []dx;
    delete []coeff_reac;
    delete []coeff_prod;
    delete []coeff_totl;
    delete []coeff_actv;
}

void
ReactionDiffusion::_fill_local_r(const double * const restrict yi,
				 double * const restrict local_r) const
{
    // intent(out) :: local_r
    for (int rxni=0; rxni<nr; ++rxni){
        // reaction rxni
        local_r[rxni] = k[rxni];
        for (int rnti=0; rnti<stoich_actv[rxni].size(); ++rnti){
            // reactant index rnti
            int si = stoich_actv[rxni][rnti];
            local_r[rxni] *= yi[si];
        }
    }
}

// The indices of x, fluxes and bins
// <indices.png>


#define C(i) y[(i)*n+si]
double
ReactionDiffusion::flux(int bi, int si, const double * const restrict y) const
{
    // bi: bin index, si: species index
    return -D[si]*(C(bi+1) - C(bi))/dx[bi];
}
#undef C

// 4π cancel
#define L(i) (x[i+1]-x[i])
#define A(i) (x[i]*x[i]) // 4πrₖ²
#define V(i) ((x[i+1]*x[i+1]*x[i+1] - x[i]*x[i]*x[i])/3) // 4πrₖ₊₁³/3 - 4πrₖ³/3
//define FLUX(i) fluxes[(i)+si*n]
#define FLUX(i) fluxes[(i)*n+si]
double
ReactionDiffusion::diffusion_contrib(int bi, int si, const double * const restrict fluxes) const
{
    // bi: bin index, si: species index, fluxes: mol/m2/s through right wall of bin
    double contrib = 0;
    switch(geom){
    case Geom::FLAT :
        if (bi > 0)   contrib += FLUX(bi-1)/L(bi);
        if (bi < N-1) contrib -= FLUX(bi)/L(bi);
        break;
    case Geom::SPHERICAL :
        if (bi > 0)   contrib += FLUX(bi-1)*A(bi)/V(bi);
        if (bi < N-1) contrib -= FLUX(bi)*A(bi+1)/V(bi);
        break;
    }
    return contrib;
}

double
ReactionDiffusion::diffusion_contrib_jac_prev(int bi) const
{
    switch(geom){
    case Geom::FLAT :      return 1.0/dx[bi-1]/L(bi);
    case Geom::SPHERICAL : return 1.0/dx[bi-1]*A(bi)/V(bi);
    }
    return 0.0/0.0; // Cylindrical goes here later
}

double
ReactionDiffusion::diffusion_contrib_jac_next(int bi) const
{
    switch(geom){
    case Geom::FLAT :      return 1.0/dx[bi]/L(bi);
    case Geom::SPHERICAL : return 1.0/dx[bi]*A(bi+1)/V(bi);
    }
    return 0.0/0.0; // Cylindrical goes here later
}
#undef L
#undef A
#undef V


#define DCDT dydt[bi*n+si]
void
ReactionDiffusion::f(double t, const double * const restrict y, double * const restrict dydt) const
{
    // initialize fluxes
    double * const fluxes = new double[(N-1)*n];
    ${"#pragma omp parallel for" if USE_OPENMP else ""}
    for (int bi=0; bi<N-1; ++bi)
        for (int si=0; si<n; ++si)
            FLUX(bi) = flux(bi, si, y);

    ${"double * local_r = new double[nr];" if not USE_OPENMP else ""}
    ${"#pragma omp parallel for" if USE_OPENMP else ""}
    for (int bi=0; bi<N; ++bi){
        // compartment bi
        ${"double * local_r = new double[nr];" if USE_OPENMP else ""}

        for (int si=0; si<n; ++si)
            DCDT = 0.0; // zero out

        // Contributions from reactions
        // ----------------------------
	const double * const yi = y+bi*n;
        _fill_local_r(yi, local_r);
        for (int rxni=0; rxni<nr; ++rxni){
            // reaction index rxni
            for (int si=0; si<n; ++si){
                // species index si
                int overall = coeff_totl[rxni*n + si];
                if (overall != 0)
                    DCDT += overall*local_r[rxni];
            }
        }

        if (N>1){
            // Contributions from diffusion
            // ----------------------------
#ifdef DEBUG2
	    printf("Fluxes:\n");
	    for (int i=0; i<N-1; ++i)
		printf("%12.5e ",fluxes[i]);
	    printf("\n");
#endif
            for (int si=0; si<n; ++si)// species index si
#ifdef DEBUG2
	    {
		double tmp = diffusion_contrib(bi, si, fluxes);
		DCDT+=tmp;
		printf("diffusion_contrib(%d, %d, fluxes)=%12.5e\n",bi,si,tmp);
	}
#else
                DCDT += diffusion_contrib(bi, si, fluxes);
#endif
        }

        ${"delete []local_r;" if USE_OPENMP else ""}

    }

    ${"delete []local_r;" if not USE_OPENMP else ""}
    delete []fluxes;
}
#undef DCDT
#undef FLUX



%for token, imaj, imin in						\
    [('dense_jac_rmaj', '(bri)*n + ri', '(bci)*n + ci'),		\
     ('dense_jac_cmaj', '(bci)*n + ci', '(bri)*n + ri'),		\
     ('banded_jac_cmaj','(bci)*n + ci', '2*n+(bri)*n+ri-((bci)*n+ci)'),	\
     ('banded_packed_jac_cmaj','(bci)*n+ci','n+(bri)*n+ri-(bci)*n-ci'), \
        ]:
#define JAC(bri, bci, ri, ci) ja[(${imaj})*ldj+${imin}]
void
ReactionDiffusion::${token}(double t, const double * const restrict y,
                            double * const restrict ja, int ldj) const
{
    ${"double * local_r = new double[nr];" if not USE_OPENMP else ""}
    ${"#pragma omp parallel for" if USE_OPENMP else ""}
    for (int bi=0; bi<N; ++bi){
        const double * const C = y + bi*n; // Conc. in `bi:th` compartment

	${'double * local_r = new double[nr];' if USE_OPENMP else ''}
    
        // Contributions from reactions
        // ----------------------------
        _fill_local_r(C, local_r);
        for (int si=0; si<n; ++si){
            // species si
            for (int dsi=0; dsi<n; ++dsi){
                // derivative wrt species dsi
                // j_i[si, dsi] = Sum_l(n_lj*Derivative(r[l], C[dsi]))
                JAC(bi,bi,si,dsi) = 0.0;
                for (int rxni=0; rxni<nr; ++rxni){
                    // reaction rxni
                    if (coeff_totl[rxni*n + si] == 0)
                        continue; // si unaffected by reaction
                    else if (coeff_actv[rxni*n + dsi] == 0)
                        continue; // rate of reaction unaffected by dsi
                    JAC(bi,bi,si,dsi) += coeff_totl[rxni*n + si]*\
                        coeff_actv[rxni*n + dsi]*local_r[rxni]/C[dsi];
                }
            }
        }

        if (N>1){
            // Contributions from diffusion
            // ----------------------------
	    if (bi > 0){
		double tmp = diffusion_contrib_jac_prev(bi);
		for (int si=0; si<n; ++si){
		    // species index si
                    JAC(bi, bi-1, si, si)  = D[si]*tmp;
                    JAC(bi, bi,   si, si) -= D[si]*tmp; // from symmetry
                }
	    }
	    if (bi < N-1){
		double tmp = diffusion_contrib_jac_next(bi);
		for (int si=0; si<n; ++si){
		    // species index si
                    JAC(bi, bi+1, si, si)  = D[si]*tmp;
                    JAC(bi, bi,   si, si) -= D[si]*tmp; // from symmetry
                }
            }
        }
	${'delete []local_r;' if USE_OPENMP else ''}
    }

    ${'delete []local_r;' if not USE_OPENMP else ''}
}
#undef JAC
%endfor
