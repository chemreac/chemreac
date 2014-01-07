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
%else:
#define omp_get_thread_num() 0
#ifdef _OPENMP
#include "CRASH_COMPILER_PLEASE_OPENMP_NOT_ENABLED"
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
    int geom
    ):
    n(n), N(N), stoich_reac(stoich_reac), 
    stoich_prod(stoich_prod), stoich_actv(stoich_actv),
    k(k), D(D), x(x), mode(mode), geom(geom)
{
    if (stoich_reac.size() != stoich_prod.size())
	throw std::length_error(
	    "stoich_reac and stoich_prod of different sizes.");
    if (D.size() != n)
	throw std::length_error(
	    "Number of diff. coeff. does not match number of species.");
    if (x.size() != N + 1)
	throw std::length_error(
	    "Number bin edges != number of compartments + 1.");

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
    delete []coeff_reac;
    delete []coeff_prod;
    delete []coeff_totl;
}

void
ReactionDiffusion::_fill_local_r(double * yi, double * local_r)
{
    // intent(out) :: local_r
    for (int j=0; j<nr; ++j){
	// reaction j
	local_r[j] = k[j];
	for (int l=0; l<stoich_reac[j].size(); ++l){
	    // reactant index l
	    int m = stoich_reac[j][l];
	    local_r[j] *= yi[m];
	}
    }
}

double *
ReactionDiffusion::_get_p(int i, double * y)
{
    // Concentrations in previous compartments (mind the boundary)
    if (i == 0)
	return y + i*n; // No diffusion to the left (solid interface)
    else
	return y + (i-1)*n;
}

double *
ReactionDiffusion::_get_n(int i, double * y)
{
    // Concentrations in previous compartments (mind the boundary)
    if (i == N-1)
	return y + i*n; // No diffusion to the right (solid interface)
    else
	return y + (i+1)*n; 
}

void
ReactionDiffusion::_get_dx(int i, double * dx_p, 
			   double * dx, double * dx_n)
{
    if (i == 0){
	*dx_n = x[i+1]-x[i];
	*dx_p = *dx_n;
    }else if (i == N-1){
	*dx_p = x[i]-x[i-1];
	*dx_n = *dx_p;
    }else{
	*dx_p = x[i]-x[i-1];
	*dx_n = x[i+1]-x[i];
    }
    *dx = (*dx_p + *dx_n)/2.0;
}

void
ReactionDiffusion::f(double t, double * y, double * dydt)
{
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
    
	// Contributions from diffusion
	// ----------------------------
	double * p_ = _get_p(i, y);
	double * n_ = _get_n(i, y);
	double * C = y + i*n;
	double diffusion_contrib = 0.0;
	double dx_p, dx, dx_n;
	_get_dx(i, &dx_p, &dx, &dx_n);
#ifdef DEBUG
	printf("dx_p=%12.5e\n", dx_p);
	printf("dx=%12.5e\n", dx);
	printf("dx_n=%12.5e\n", dx_n);
#endif
	for (int k=0; k<n; ++k){
	    if (i > 0)
		diffusion_contrib += SA_p*(p_[k] - C[k])/dx_p/V;
	    if (i < N-1)
		diffusion_contrib += SA_n*(n_[k] - C[k])/dx_n/V;
	    diffusion_contrib *= D[k];
	    dydt[i*n + k] += diffusion_contrib;
	}
	${"delete []local_r;" if USE_OPENMP else ""}

    }

    ${"delete []local_r;" if not USE_OPENMP else ""}
    nfeval++;
}


%for token, imaj, imin in					\
    [('dense_jac_rmaj', '((bri)*n + ri)', '(bci)*n + ci'), \
     ('dense_jac_cmaj', '((bci)*n + ci)', '(bri)*n + ri'),\
     ('banded_jac_cmaj','((bci)*n + ci)', '2*n+((bri)*n+ri)-((bci)*n+ci)'),\
     ('banded_packed_jac_cmaj','((bci)*n + ci)', 'n+((bri)*n+ri)-((bci)*n+ci)'),\
	]:
#define JAC(bri, bci, ri, ci) ja[${imaj}*ldj+${imin}]
    void
    ReactionDiffusion::${token}(double t, double * y, double * ja, 
				int ldj)
    {
	// intent(out) :: ja
	// Assume ja zeroed out on entry

#ifndef _OPENMP
	double * local_r = new double[nr];
#else
#pragma omp parallel for
#endif
	for (int i=0; i<N; ++i){
	    double * C = y + i*n; // Conc. in `i:th` compartment

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
		    JAC(i,i,j,m) = 0.0;
		    for (int l=0; l<nr; ++l){
			// reaction l
			if (coeff_totl[l*n + j] == 0)
			    continue;
			else if (coeff_reac[l*n + m] == 0)
			    continue;
			JAC(i,i,j,m) += coeff_totl[l*n + j]*\
			    coeff_reac[l*n + m]*local_r[l]/C[m];
		    }
		}
	    }

	    double dx_p, dx, dx_n;
	    _get_dx(i, &dx_p, &dx, &dx_n);
	    // Contributions from diffusion
	    // ----------------------------
	    for (int j=0; j<n; ++j){
		// species j
		if (i > 0){ // Subdiagonal (diffusion from left)
		    JAC(i,i-1,j,j)  = SA_p*D[j]/dx_p/V; 
		    JAC(i,  i,j,j) -= SA_p*D[j]/dx_p/V;
		}
		if (i < N-1){ // Superdiagonal (diffusion from right)
		    JAC(i,i+1,j,j)  = SA_n*D[j]/dx_n/V; 
		    JAC(i,  i,j,j) -= SA_n*D[j]/dx_n/V;
		}
	    }
#ifdef _OPENMP
	    delete []local_r;
#endif
	}

#ifndef _OPENMP
	delete []local_r;
#endif
	njeval++;
    }
#undef JAC
%endfor
