## -*- coding: utf-8 -*-
// ${_warning_in_the_generated_file_not_to_edit}
<%doc>
// This is a templated source file.
// Render template using Mako (Python templating engine)
</%doc>
#include <algorithm> // count
#include <vector>
#include "chemreac.h"

#ifdef DEBUG
#include <cstdio>
#include <iostream>
#define PRINT_ARR(ARR, LARR) for(int i_=0; i_<LARR; ++i_) {std::cout << ARR[i_] << " " << std::endl;};
#endif

%if USE_OPENMP:
#ifndef _OPENMP
  #error
#endif
#include <omp.h>
%else:
#ifdef _OPENMP
  #error
#endif
#define omp_get_thread_num() 0
%endif

#ifndef NULL
#define NULL 0
#endif // NULL not part of C++ standard

using std::vector;
using std::count;
namespace chemreac {

#define D_JAC(bi, j) D_jac[3*(bi) + j]
#define D_WEIGHT(bi, j) D_weight[2*(bi) + j]
// 1D discretized reaction diffusion
ReactionDiffusion::ReactionDiffusion(
    int n,
    vector<vector<int> > stoich_reac,
    vector<vector<int> > stoich_prod,
    vector<double> k,
    int N, 
    vector<double> D,
    vector<double> x, // separation
    vector<vector<int> > stoich_actv_, // vectors of size 0 in stoich_actv_ => "copy from stoich_reac"
    vector<vector<double> > bin_k_factor, // per bin modulation of first k's
    vector<int> bin_k_factor_span, // modulation over reactions
    int geom_,
    bool logy,
    bool logt
    ):
    n(n), stoich_reac(stoich_reac), stoich_prod(stoich_prod),
    k(k), N(N), D(D), x(x), bin_k_factor(bin_k_factor), 
    bin_k_factor_span(bin_k_factor_span), logy(logy), logt(logt), nr(stoich_reac.size())
{
    if (N == 2) throw std::logic_error("2nd order PDE requires at least 3 stencil points.");
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

    D_weight = new double[2*N];
    D_jac = new double[3*N];
    switch(geom_) { // Precalc coeffs for Jacobian for current geom.
  %for IDX, GEOMSTR in enumerate(['FLAT', 'CYLINDRICAL', 'SPHERICAL']):
    case ${IDX}:
        geom = Geom::${GEOMSTR};
        if (N == 1) break;
        {
            int i = 1;
            D_WEIGHT(0,0) = ${COEFF[GEOMSTR,-1,'bw']};
            D_WEIGHT(0,1) = ${COEFF[GEOMSTR,-1,'fw']};
            D_JAC(0,1) = ${COEFF[GEOMSTR,-1,'Jc']};
            D_JAC(0,2) = ${COEFF[GEOMSTR,-1,'Jn']};           
            do {
                D_WEIGHT(i,0) = ${COEFF[GEOMSTR,0,'bw']};
                D_WEIGHT(i,1) = ${COEFF[GEOMSTR,0,'fw']};
                D_JAC(i,0) = ${COEFF[GEOMSTR,0,'Jp']};
                D_JAC(i,1) = ${COEFF[GEOMSTR,0,'Jc']};           
                D_JAC(i,2) = ${COEFF[GEOMSTR,0,'Jn']};           
            } while (i++ < N-2);
            i--;
            #ifdef DEBUG
            printf("i=%d, N=%d\n", i, N);
            #endif
            D_WEIGHT(i+1,0) = ${COEFF[GEOMSTR,1,'bw']};
            D_WEIGHT(i+1,1) = ${COEFF[GEOMSTR,1,'fw']};            
            D_JAC(i+1, 0) = ${COEFF[GEOMSTR,1,'Jp']};
            D_JAC(i+1, 1) = ${COEFF[GEOMSTR,1,'Jc']};           
        }
        break;
  %endfor
    default: throw std::logic_error("Unknown geom.");
    }

    //    nr = stoich_reac.size();

    for (int ri=0; ri<nr; ++ri){
        for (vector<int>::iterator si=stoich_reac[ri].begin(); si != stoich_reac[ri].end(); ++si)
            if (*si > n-1)
                throw std::logic_error("At least one species index in stoich_reac > (n-1)");
        for (vector<int>::iterator si=stoich_prod[ri].begin(); si != stoich_prod[ri].end(); ++si)
            if (*si > n-1)
                throw std::logic_error("At least one species index in stoich_prod > (n-1)");
        for (vector<int>::iterator si=stoich_actv_[ri].begin(); si != stoich_actv_[ri].end(); ++si)
            if (*si > n-1)
                throw std::logic_error("At least one species index in stoich_actv > (n-1)");
    }

    coeff_reac = new int[nr*n];
    coeff_prod = new int[nr*n];
    coeff_totl = new int[nr*n];
    coeff_actv = new int[nr*n];

    stoich_actv.reserve(nr);
    for (int rxni=0; rxni<nr; ++rxni){ // reaction index 
        if (stoich_actv_[rxni].size() == 0)
            stoich_actv.push_back(stoich_reac[rxni]); // massaction
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

    // Handle bin_k_factors:
    for (int i=0; i<bin_k_factor_span.size(); ++i)
        for (int j=0; j<bin_k_factor_span[i]; ++j)
            i_bin_k.push_back(i);
    n_factor_affected_k = i_bin_k.size();
}

ReactionDiffusion::~ReactionDiffusion()
{
    delete []D_weight;
    delete []D_jac;
    delete []coeff_reac;
    delete []coeff_prod;
    delete []coeff_totl;
    delete []coeff_actv;
}

#define FACTOR(ri, bi) (((ri) < n_factor_affected_k) ? \
            bin_k_factor[bi][i_bin_k[ri]] : 1)
void
ReactionDiffusion::_fill_local_r(int bi, const double * const restrict y,
                 double * const restrict local_r) const
{
    // intent(out) :: local_r
    for (int rxni=0; rxni<nr; ++rxni){
        // reaction rxni
        if (logy)
            local_r[rxni] = 0;
        else
            local_r[rxni] = 1;

        for (int rnti=0; rnti<stoich_actv[rxni].size(); ++rnti){
            // reactant index rnti
            int si = stoich_actv[rxni][rnti];
            if (logy)
                local_r[rxni] += y[bi*n+si];
            else
                local_r[rxni] *= y[bi*n+si];
        }

        if (logy)
            local_r[rxni] = exp(local_r[rxni]);

        local_r[rxni] *= FACTOR(rxni,bi)*k[rxni];
    }
}
#undef FACTOR

// The indices of x, fluxes and bins
// <indices.png>

#define Y(bi, si) y[(bi)*n+(si)]
#define C(bi, si) ( (logy) ? exp(Y(bi, si)) : Y(bi, si) )
#define DCB(bi, si) ((bi>0) ? (C(bi,si) - C(bi-1,si)) : 0)
#define DCF(bi, si) ((bi<N-1) ? (C(bi,si+1) - C(bi,si)) : 0)
#define DYDT(bi, si) dydt[(bi)*(n)+(si)]
void
ReactionDiffusion::f(double t, const double * const restrict y, double * const restrict dydt) const
{
    ${"double * local_r = new double[nr];" if not USE_OPENMP else ""}
    ${"#pragma omp parallel for" if USE_OPENMP else ""}
    for (int bi=0; bi<N; ++bi){
        // compartment bi
        ${"double * local_r = new double[nr];" if USE_OPENMP else ""}

        for (int si=0; si<n; ++si)
            DYDT(bi, si) = 0.0; // zero out

        // Contributions from reactions
        // ----------------------------
        _fill_local_r(bi, y, local_r);
        for (int rxni=0; rxni<nr; ++rxni){
            // reaction index rxni
            for (int si=0; si<n; ++si){
                // species index si
                int overall = coeff_totl[rxni*n + si];
                if (overall != 0)
                    DYDT(bi, si) += overall*local_r[rxni];
            }
        }
        if (N>1){
            // Contributions from diffusion
            // ----------------------------
            for (int si=0; si<n; ++si){ // species index si
                if (D[si] != 0.0)
                    DYDT(bi, si) += D[si]*(D_WEIGHT(bi, 0)*DCB(bi, si)+\
                                           D_WEIGHT(bi, 1)*DCF(bi, si));
            }
        }
        if (logy){
            if (logt)
                for (int si=0; si<n; ++si)
                    DYDT(bi, si) *= exp(t-Y(bi,si));
            else
                for (int si=0; si<n; ++si)
                    DYDT(bi, si) *= exp(-Y(bi,si));
        } else {
            if (logt)
                for (int si=0; si<n; ++si)
                    DYDT(bi, si) *= exp(t);
        }

        ${"delete []local_r;" if USE_OPENMP else ""}

    }

    ${"delete []local_r;" if not USE_OPENMP else ""}
}
#undef DCB
#undef DCF
#undef DCDT // Y(bi, si) still defined.
#undef D_WEIGHT


%for token, imaj, imin in [\
    ('dense_jac_rmaj',         '(bri)*n+ri', '(bci)*n + ci'),\
    ('dense_jac_cmaj',         '(bci)*n+ci', '(bri)*n + ri'),\
    ('banded_packed_jac_cmaj', '(bci)*n+ci', '(1+bri-(bci))*n+ri-ci'),\
    ('banded_padded_jac_cmaj', '(bci)*n+ci', '(2+bri-(bci))*n+ri-ci'),\
    ]:
#define JAC(bri, bci, ri, ci) ja[(${imaj})*ldj+${imin}]
void
ReactionDiffusion::${token}(double t, const double * const restrict y,
                            double * const restrict ja, int ldj) const
{
    // `t`: time (log(t) if logt=1)
    // `y`: concentrations (log(conc) if logy=True)
    // `ja`: jacobian (allocated 1D array to hold dense or banded)
    // `ldj`: leading dimension of ja (useful for padding)
    double * restrict fout = NULL;
    if (logy){
        fout = new double[n*N];
        f(t, y, fout);
    }

    ${'double * local_r = new double[nr];' if not USE_OPENMP else ''}
    ${'#pragma omp parallel for' if USE_OPENMP else ''}
    for (int bi=0; bi<N; ++bi){
        // Conc. in `bi:th` compartment
        ${'double * local_r = new double[nr];' if USE_OPENMP else ''}
    
        // Contributions from reactions
        // ----------------------------
        _fill_local_r(bi, y, local_r);
        for (int si=0; si<n; ++si){
            // species si
            for (int dsi=0; dsi<n; ++dsi){
                // derivative wrt species dsi
                // j_i[si, dsi] = Sum_l(n_lj*Derivative(r[l], local_y[dsi]))
                JAC(bi,bi,si,dsi) = 0.0;
                for (int rxni=0; rxni<nr; ++rxni){
                    // reaction rxni
                    if (coeff_totl[rxni*n + si] == 0)
                        continue; // si unaffected by reaction
                    if (coeff_actv[rxni*n + dsi] == 0)
                        continue; // rate of reaction unaffected by dsi
                    double tmp = coeff_totl[rxni*n + si]*\
                    coeff_actv[rxni*n + dsi]*local_r[rxni];
                    if (!logy)
                        tmp /= Y(bi,dsi);
                    JAC(bi,bi,si,dsi) += tmp;
                }
                if (logy)
                    JAC(bi,bi,si,dsi) *= exp(-Y(bi,si));
            }
        }

        // Contributions from diffusion
        // ----------------------------
        if (bi > 0){
            // Diffusion over left boundary
            for (int si=0; si<n; ++si){
                // species index si
                if (D[si] == 0.0) continue;
                JAC(bi, bi, si, si) += D[si]*D_JAC(bi, 1);
                if (logy)
                    JAC(bi, bi-1, si, si)  = D[si]*D_JAC(bi, 0)*exp(Y(bi-1,si)-Y(bi,si));
                else
                    JAC(bi, bi-1, si, si)  = D[si]*D_JAC(bi, 0);
            }
        }
        if (bi < N-1){
            // Diffusion over right boundary
            for (int si=0; si<n; ++si){
                // species index si
                if (D[si] == 0.0) continue;
                JAC(bi, bi, si, si) += D[si]*D_JAC(bi, 1);
                if (logy)
                    JAC(bi, bi+1, si, si)  = D[si]*D_JAC(bi, 2)*exp(Y(bi+1,si)-Y(bi,si));
                else
                    JAC(bi, bi+1, si, si)  = D[si]*D_JAC(bi, 2);
            }
        }

        // Logartihmic time
        // ----------------------------
        if (logt){
            for (int si=0; si<n; ++si){
                for (int dsi=0; dsi<n; ++dsi){
                    JAC(bi, bi, si, dsi) *= exp(t);
                }
                if (bi>0)
                    JAC(bi, bi-1, si, si) *= exp(t);
                if (bi<N-1)
                    JAC(bi, bi+1, si, si) *= exp(t);
            }
        }

        // Logrithmic concentrations
        // ----------------------------
        if (logy){
            for (int si=0; si<n; ++si)
                JAC(bi, bi, si, si) -= fout[bi*n+si];
        }

        ${'delete []local_r;' if USE_OPENMP else ''}
    }
    ${'delete []local_r;' if not USE_OPENMP else ''}
    if (logy)
        delete []fout;
}
#undef JAC
%endfor
#undef D_JAC
#undef Y

void ReactionDiffusion::per_rxn_contrib_to_fi(double t, const double * const restrict y,
                                              int si, double * const restrict out) const
{
    double * local_r = new double[nr];
    _fill_local_r(0, y, local_r);
    for (int ri=0; ri<nr; ++ri){
	out[ri] = coeff_totl[ri*n+si]*local_r[ri];
    }
    delete []local_r;
}

int ReactionDiffusion::get_geom_as_int() const
{
    switch(geom){
    case Geom::FLAT :        return 0;
    case Geom::CYLINDRICAL : return 1;
    case Geom::SPHERICAL :   return 2;
    default:                 return -1;
    }
}

}; // namespace chemreac
