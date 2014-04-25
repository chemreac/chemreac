## -*- coding: utf-8 -*-
// ${_warning_in_the_generated_file_not_to_edit}
<%doc>
// This is a templated source file.
// Render template using Mako (Python templating engine)
</%doc>
#include <algorithm> // std::count
#include <vector>    // std::vector
#include <algorithm> // std::max, std::min
#include <cstdlib> // free,  C++11 aligned_alloc
#include "chemreac.h"
#include "c_fornberg.h" // fintie differences

#ifdef DEBUG
#include <cstdio>
#include <iostream>
#define PRINT_ARR(ARR, LARR) for(int i_=0; i_<LARR; ++i_) {std::cout << ARR[i_] << " ";}; std::cout << std::endl;
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

using std::vector;
using std::count;
using std::min;
using std::max;

namespace chemreac {

// 1D discretized reaction diffusion
ReactionDiffusion::ReactionDiffusion(
    uint n,
    const vector<vector<uint> > stoich_reac,
    const vector<vector<uint> > stoich_prod,
    vector<double> k,
    uint N, 
    vector<double> D,
    const vector<double> x, // separation
    vector<vector<uint> > stoich_actv_, // vectors of size 0 in stoich_actv_ => "copy from stoich_reac"
    vector<vector<double> > bin_k_factor, // per bin modulation of first k's
    vector<uint> bin_k_factor_span, // modulation over reactions
    int geom_,
    bool logy,
    bool logt,
    uint nstencil,
    bool lrefl,
    bool rrefl
    ):
    n(n), N(N), nstencil(nstencil), nr(stoich_reac.size()),
    logy(logy), logt(logt), stoich_reac(stoich_reac), stoich_prod(stoich_prod),
    k(k),  D(D), x(x), bin_k_factor(bin_k_factor), 
    bin_k_factor_span(bin_k_factor_span), lrefl(lrefl), rrefl(rrefl)
{
    if (N == 0) throw std::logic_error("Zero bins sounds boring.");
    if (N == 2) throw std::logic_error("2nd order PDE requires at least 3 stencil points.");
    if (nstencil % 2 == 0) throw std::logic_error("Only odd number of stencil points supported");
    if ((N == 1) && (nstencil != 1)) throw std::logic_error("You must set nstencil=1 for N=1");
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
    case 0:
        geom = Geom::FLAT;
        break;
    case 1:
        geom = Geom::CYLINDRICAL;
        break;
    case 2:
        geom = Geom::SPHERICAL;
        break;
    default:
        throw std::logic_error("Unknown geom.");
    }

    // Finite difference scheme
    xc = new double[N+nstencil-1];
    for (uint i=0; i<N; ++i) xc[i+(nstencil-1)/2] = (x[i+1] + x[i])/2;
    for (uint i=0; i<(nstencil-1)/2; ++i){
        // reflection
        int nsidep = (nstencil-1)/2;
        xc[nsidep-i-1] = xc[nsidep-i] - (xc[nsidep-i+1] - xc[nsidep-i]);
        xc[N+nsidep+i] = xc[N+nsidep+i-1] + (xc[N+nsidep+i-1] - xc[N+nsidep+i-2]);
    }
    D_weight = new double[nstencil*N];

    // Precalc coeffs for Jacobian for current geom.
    // not centered diffs close to boundaries
    for (uint bi=0; bi<N; bi++)
        _apply_fd(bi);

    for (uint ri=0; ri<nr; ++ri){
        for (auto si=stoich_reac[ri].begin(); si != stoich_reac[ri].end(); ++si)
            if (*si > n-1)
                throw std::logic_error("At least one species index in stoich_reac > (n-1)");
        for (auto si=stoich_prod[ri].begin(); si != stoich_prod[ri].end(); ++si)
            if (*si > n-1)
                throw std::logic_error("At least one species index in stoich_prod > (n-1)");
        for (auto si=stoich_actv_[ri].begin(); si != stoich_actv_[ri].end(); ++si)
            if (*si > n-1)
                throw std::logic_error("At least one species index in stoich_actv > (n-1)");
    }

    coeff_reac = new int[nr*n];
    coeff_prod = new int[nr*n];
    coeff_totl = new int[nr*n];
    coeff_actv = new int[nr*n];

    stoich_actv.reserve(nr);
    for (uint rxni=0; rxni<nr; ++rxni){ // reaction index 
        if (stoich_actv_[rxni].size() == 0)
            stoich_actv.push_back(stoich_reac[rxni]); // massaction
    else
        stoich_actv.push_back(stoich_actv_[rxni]);
        for (uint si=0; si<n; ++si){ // species index
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
    for (uint i=0; i<bin_k_factor_span.size(); ++i)
        for (uint j=0; j<bin_k_factor_span[i]; ++j)
            i_bin_k.push_back(i);
    n_factor_affected_k = i_bin_k.size();
}

ReactionDiffusion::~ReactionDiffusion()
{
    delete []xc;
    delete []D_weight;
    delete []coeff_reac;
    delete []coeff_prod;
    delete []coeff_totl;
    delete []coeff_actv;
}


#define D_WEIGHT(bi, j) D_weight[nstencil*(bi) + j]
#define FDWEIGHT(order, local_index) c[nstencil*(order) + local_index]
void ReactionDiffusion::_apply_fd(uint bi){
    double * c = new double[3*nstencil];
    double * lxc = new double[nstencil]; // local shifted x-centers
    uint nsidep = (nstencil-1)/2;
    uint around = bi + nsidep;
    uint start  = bi;
    if (!lrefl) // shifted finite diff
        start = max(nsidep, start);
    if (!rrefl) // shifted finite diff
        start = min(N-nstencil+nsidep, bi);

    for (uint li=0; li<nstencil; ++li) // li: local index
        lxc[li] = xc[start + li] - xc[around];
    fornberg_populate_weights(0, lxc, nstencil-1, 2, c);
    delete []lxc;

    for (uint li=0; li<nstencil; ++li){ // li: local index
        D_WEIGHT(bi, li) = FDWEIGHT(2, li);
        switch(geom){
        case Geom::CYLINDRICAL: // Laplace operator in cyl coords.
            D_WEIGHT(bi, li) += FDWEIGHT(1, li)*1/xc[around];
            break;
        case Geom::SPHERICAL: // Laplace operator in sph coords.
            D_WEIGHT(bi, li) += FDWEIGHT(1, li)*2/xc[around];
            break;
        default:
            break;
        }
    }
    delete []c;
}
#undef FDWEIGHT

#define FACTOR(ri, bi) ( ((ri) < n_factor_affected_k) ? \
            bin_k_factor[bi][i_bin_k[ri]] : 1 )
void
ReactionDiffusion::_fill_local_r(int bi, const double * const restrict y,
                 double * const restrict local_r) const
{
    // intent(out) :: local_r
    for (uint rxni=0; rxni<nr; ++rxni){
        // reaction rxni
        if (logy)
            local_r[rxni] = 0;
        else
            local_r[rxni] = 1;

        for (uint rnti=0; rnti<stoich_actv[rxni].size(); ++rnti){
            // reactant index rnti
            int si = stoich_actv[rxni][rnti];
            if (logy)
                local_r[rxni] += y[bi*n+si];
            else
                local_r[rxni] *= y[bi*n+si];
        }

        local_r[rxni] *= FACTOR(rxni,bi)*k[rxni];
        if (logy)
            local_r[rxni] = exp(local_r[rxni]);
    }
}
#undef FACTOR

// The indices of x, fluxes and bins
// <indices.png>

#define Y(bi, si) y[(bi)*n+(si)]
#define LC(bi, si) liny[(bi)*n+(si)]
#define DYDT(bi, si) dydt[(bi)*(n)+(si)]
void
ReactionDiffusion::f(double t, const double * const restrict y, double * const restrict dydt) const
{
    double * liny __attribute__((aligned(16))) = nullptr;
    if ((N > 1) && logy){
        int nliny = 16/sizeof(double)*int(ceil(sizeof(double)*n*N / 16.0));
        liny = (double * const)(aligned_alloc(16,nliny*sizeof(double)));

        ${"#pragma omp parallel for if (N > 2)" if USE_OPENMP else ""}
        for (uint bi=0; bi<N; ++bi)
            for (uint si=0; si<n; ++si)
                LC(bi, si) = exp(Y(bi, si));
    }

    ${"double * local_r = new double[nr];" if not USE_OPENMP else ""}
    ${"#pragma omp parallel for if (N > 2)" if USE_OPENMP else ""}
    for (uint bi=0; bi<N; ++bi){
        // compartment bi
        ${"double * local_r = new double[nr];" if USE_OPENMP else ""}

        for (uint si=0; si<n; ++si)
            DYDT(bi, si) = 0.0; // zero out

        // Contributions from reactions
        // ----------------------------
        _fill_local_r(bi, y, local_r);
        for (uint rxni=0; rxni<nr; ++rxni){
            // reaction index rxni
            for (uint si=0; si<n; ++si){
                // species index si
                int overall = coeff_totl[rxni*n + si];
                if (overall != 0)
                    DYDT(bi, si) += overall*local_r[rxni];
            }
        }
        if (N>1){
            // Contributions from diffusion
            // ----------------------------
            uint nsidep = (nstencil-1)/2;
            int starti;
            if ((bi < nsidep) && (!lrefl)){
                starti = 0;
            } else if ((bi >= N-nsidep) && (!rrefl)){
                starti = N-nstencil;
            } else{
                starti = bi-nsidep;
            }
            for (uint si=0; si<n; ++si){ // species index si
                if (D[si] == 0.0) continue;
                double tmp = 0;
                for (uint xi=0; xi<nstencil; ++xi){
                    int biw = starti + xi;
                    // reflective logic:
                    if (starti < 0)
                        biw = abs(biw); // lrefl==true
                    else if (starti >= N-nsidep)
                        biw = N - 1 - abs(N-1-biw); // rrefl==true

                    tmp += D_WEIGHT(bi, xi) * ((logy) ? LC(biw, si) : Y(biw, si));
                }
                DYDT(bi, si) += D[si]*tmp;
            }
        }
        if (logy){
            if (logt)
                for (uint si=0; si<n; ++si)
                    DYDT(bi, si) *= exp(t-Y(bi,si));
            else
                for (uint si=0; si<n; ++si)
                    DYDT(bi, si) *= exp(-Y(bi,si));
        } else {
            if (logt)
                for (uint si=0; si<n; ++si)
                    DYDT(bi, si) *= exp(t);
        }

        ${"delete []local_r;" if USE_OPENMP else ""}

    }
    ${"delete []local_r;" if not USE_OPENMP else ""}

    if ((N > 2) && logy)
        free((void*)liny);
}
#undef DCDT 
// Y(bi, si) and LC(bi, si) still defined.


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
    // Note: does not return a strictly correct Jacobian, only 1 pair of bands.
    // `t`: time (log(t) if logt=1)
    // `y`: concentrations (log(conc) if logy=True)
    // `ja`: jacobian (allocated 1D array to hold dense or banded)
    // `ldj`: leading dimension of ja (useful for padding)
    double * restrict fout = nullptr;
    if (logy){
        fout = new double[n*N];
        f(t, y, fout);
    }

    ${'double * local_r = new double[nr];' if not USE_OPENMP else ''}
    ${'#pragma omp parallel for' if USE_OPENMP else ''}
    for (uint bi=0; bi<N; ++bi){
        // Conc. in `bi:th` compartment
        ${'double * local_r = new double[nr];' if USE_OPENMP else ''}
    
        // Contributions from reactions
        // ----------------------------
        _fill_local_r(bi, y, local_r);
        for (uint si=0; si<n; ++si){
            // species si
            for (uint dsi=0; dsi<n; ++dsi){
                // derivative wrt species dsi
                // j_i[si, dsi] = Sum_l(n_lj*Derivative(r[l], local_y[dsi]))
                JAC(bi,bi,si,dsi) = 0.0;
                for (uint rxni=0; rxni<nr; ++rxni){
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
        if (N > 1) {
            // reflctive logic goes into centerli...
            int centerli = min((int)bi, max(((int)nstencil - 1)/2, (int)bi-(int)N+(int)nstencil));
            for (uint si=0; si<n; ++si){ // species index si
                if (D[si] == 0.0) continue;
                // Not a strict Jacobian only diagonal plus closest bands...
                if (!logy)
                    JAC(bi, bi, si, si) += D[si]*D_WEIGHT(bi, centerli);
                if (bi > 0)
                    JAC(bi, bi-1, si, si) = D[si]*D_WEIGHT(bi, centerli - 1)*\
                        ( (logy) ? exp(Y(bi-1, si) - Y(bi, si)) : 1 );
                if (bi < N-1)
                    JAC(bi, bi+1, si, si)  = D[si]*D_WEIGHT(bi, centerli + 1)*\
                        ( (logy) ? exp(Y(bi+1, si) - Y(bi, si)) : 1 );
            }
        }

        // Logartihmic time
        // ----------------------------
        if (logt){
            for (uint si=0; si<n; ++si){
                for (uint dsi=0; dsi<n; ++dsi){
                    JAC(bi, bi, si, dsi) *= exp(t);
                }
                if (bi>0)
                    JAC(bi, bi-1, si, si) *= exp(t);
                if (bi<N-1)
                    JAC(bi, bi+1, si, si) *= exp(t);
            }
        }
        ${'delete []local_r;' if USE_OPENMP else ''}
    }
    ${'delete []local_r;' if not USE_OPENMP else ''}
    if (logy)
        delete []fout;
}
#undef JAC
%endfor
#undef LC 
#undef Y
#undef D_WEIGHT

void ReactionDiffusion::per_rxn_contrib_to_fi(double t, const double * const restrict y,
                                              uint si, double * const restrict out) const
{
    double * local_r = new double[nr];
    _fill_local_r(0, y, local_r);
    for (uint ri=0; ri<nr; ++ri){
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
