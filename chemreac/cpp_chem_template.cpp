// This is a templated source file.
// Render template using Mako (Python templating engine)

#include <algorithm> // count
#include <vector>
#include "cpp_chem.hpp"

#ifdef _OPENMP
   #include <omp.h>
#else
   #define omp_get_thread_num() 0
#endif

using std::vector;
using std::count;

// 1D discretized reaction diffusion
ReactionDiffusion::ReactionDiffusion(int n, int N, 
				     vector<vector<int> > stoich_reac,
				     vector<vector<int> > stoich_prod,
				     vector<vector<int> > stoich_actv,
				     vector<double> k,
				     vector<double> D,
				     vector<double> x, // separation
				     int mode):
  n(n), N(N), mode(mode), stoich_reac(stoich_reac), 
  stoich_prod(stoich_prod), stoich_actv(stoich_actv), k(k), D(D), x(x)
{
  if (stoich_reac.size() != stoich_prod.size())
    throw std::length_error("stoich_reac and stoich_prod of different sizes.");
  if (D.size() != n)
    throw std::length_error("Number of diffusion coefficients does not match number of species.");
  if (x.size() != N)
    throw std::length_error("Number of separation distances does not match number of compartments.");

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
	coeff_totl[ri*n+si] = coeff_prod[ri*n+si] - coeff_reac[ri*n+si];
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
#ifdef _OPENMP
#pragma omp parallel for
#else
  double * local_r = new double[nr];
#endif
  for (int i=0; i<N; ++i){
    // compartment i
#ifdef _OPENMP
    double * local_r = new double[nr];
#endif

    for (int j=0; j<n; ++j)
      dydt[n*i+j] = 0.0; // zero out

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
    double dx_p, dx, dx_n;
    _get_dx(i, &dx_p, &dx, &dx_n);
    for (int k=0; k<n; ++k)
      dydt[i*n + k] += D[k]*(p_[k]/dx_p - 2*C[k]/dx + n_[k]/dx_n);

#ifdef _OPENMP
    delete []local_r;
#endif

  }

#ifndef _OPENMP
  delete []local_r;
#endif
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

#ifdef _OPENMP
#pragma omp parallel for
#else
  double * local_r = new double[nr];
#endif
  for (int i=0; i<N; ++i){
    double * C = y + i*n; // Concentrations in `i:th` compartment

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
	JAC(i,i-1,j,j) = D[j]/dx_p; 
	JAC(i,i,j,j) -= D[j]/dx; // Diagonal (mass balance)
      }
      if (i < N-1){ // Superdiagonal (diffusion from right)
	JAC(i,i+1,j,j) = D[j]/dx_n; 
	JAC(i,i,j,j) -= D[j]/dx; // Diagonal (mass balance)
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
