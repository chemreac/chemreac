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
#include "c_fornberg.h" // fintie differences (remember to link fortran object fornberg.o)

#ifdef DEBUG
#include <cstdio>
#include <iostream>
#define PRINT_ARR(ARR, LARR) for(int i_=0; i_<LARR; ++i_) {std::cout << ARR[i_] << " ";}; std::cout << std::endl;
#endif

%if USE_OPENMP:
#ifndef _OPENMP
  #error "Have you forgotten -fopenmp flag?"
#endif
#include <omp.h>
%else:
#ifdef _OPENMP
  #error "You should render OpenMP enabled sources"
#endif
#define omp_get_thread_num() 0
%endif

namespace chemreac {

using std::vector;
using std::count;
using std::min;
using std::max;

// 1D discretized reaction diffusion
ReactionDiffusion::ReactionDiffusion(
    uint n,
    const vector<vector<uint> > stoich_reac,
    const vector<vector<uint> > stoich_prod,
    vector<double> k,
    uint N,
    vector<double> D,
    const vector<int> z_chg,
    vector<double> mobility,
    const vector<double> x, // separation
    vector<vector<uint> > stoich_actv_, // vectors of size 0 in stoich_actv_ => "copy from stoich_reac"
    vector<vector<double> > bin_k_factor, // per bin modulation of first k's
    vector<uint> bin_k_factor_span, // modulation over reactions
    int geom_,
    bool logy,
    bool logt,
    bool logx,
    uint nstencil,
    bool lrefl,
    bool rrefl,
    bool auto_efield,
    double surf_chg,
    double eps
    ):
    n(n), N(N), nstencil(nstencil), nsidep((nstencil-1)/2), nr(stoich_reac.size()),
    logy(logy), logt(logt), logx(logx), stoich_reac(stoich_reac), stoich_prod(stoich_prod),
    k(k),  D(D), z_chg(z_chg), mobility(mobility), x(x), bin_k_factor(bin_k_factor),
    bin_k_factor_span(bin_k_factor_span), lrefl(lrefl), rrefl(rrefl), auto_efield(auto_efield),
    surf_chg(surf_chg), eps(eps), efield(new double[N])
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
        if (mobility.size() != n)
            throw std::length_error(
                "Length of mobility does not match number of species.");
        if (z_chg.size() != n)
            throw std::length_error(
                "Length of z_chg does not match number of species.");
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
    D_weight = new double[nstencil*N];
    A_weight = new double[nstencil*N];
    for (uint i=0; i<N; ++i) efield[i] = 0.0;
    xc = new double[nsidep + N + nsidep]; // xc padded with virtual bins
    for (uint i=0; i<N; ++i)
        xc[nsidep + i] = (x[i] + x[i + 1])/2;

    for (uint i=0; i<nsidep; ++i){
        // reflection
        xc[nsidep - i - 1] = 2*x[0] - xc[nsidep + i];
        xc[nsidep + i + N] = 2*x[N] - xc[nsidep + N - i - 1];
    }

    // Precalc coeffs for Jacobian for current geom.
    // not centered diffs close to boundaries
    for (uint bi=0; bi<N; bi++)
        _apply_fd(bi);

    // Stoichiometry
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
    delete []efield;
    delete []A_weight;
    delete []D_weight;
    delete []coeff_reac;
    delete []coeff_prod;
    delete []coeff_totl;
    delete []coeff_actv;
}

uint ReactionDiffusion::_stencil_bi_lbound(uint bi) const
{
    const uint le = lrefl ? 0 : nsidep;
    const uint re = rrefl ? 0 : nsidep;
    return max(le, min(N + 2*nsidep - re - nstencil, bi));
}

uint ReactionDiffusion::_xc_bi_map(uint xci) const
{
    if (xci < nsidep)
        return nsidep - xci - 1;
    else if (xci >= N+nsidep)
        return 2*N - xci;
    else
        return xci - nsidep;
}


#define D_WEIGHT(bi, li) D_weight[nstencil*(bi) + li]
#define A_WEIGHT(bi, li) A_weight[nstencil*(bi) + li]
#define FDWEIGHT(order, local_index) c[nstencil*(order) + local_index]
void ReactionDiffusion::_apply_fd(uint bi){
    double * const c = new double[3*nstencil];
    double * const lxc = new double[nstencil]; // local shifted x-centers
    uint around = bi + nsidep;
    uint start  = bi;
    if (!lrefl) // shifted finite diff
        start = max(nsidep, start);
    if (!rrefl) // shifted finite diff
        start = min(N - nstencil + nsidep, start);
    for (uint li=0; li<nstencil; ++li) // li: local index
        lxc[li] = xc[start + li] - xc[around];
    fornberg_populate_weights(0, lxc, nstencil-1, 2, c);
    delete []lxc;

    for (uint li=0; li<nstencil; ++li){ // li: local index
        D_WEIGHT(bi, li) = FDWEIGHT(2, li);
        A_WEIGHT(bi, li) = FDWEIGHT(1, li);
        if (logx){
            switch(geom){
            case Geom::FLAT:
                D_WEIGHT(bi, li) -= FDWEIGHT(1, li);
                break;
            case Geom::CYLINDRICAL:
                A_WEIGHT(bi, li) += FDWEIGHT(0, li);
                break;
            case Geom::SPHERICAL:
                D_WEIGHT(bi, li) += FDWEIGHT(1, li);
                A_WEIGHT(bi, li) += 2*FDWEIGHT(0, li);
                break;
            }
            D_WEIGHT(bi, li) *= exp(-2*xc[around]);
            A_WEIGHT(bi, li) *= exp(-xc[around]);
        } else {
            switch(geom){
            case Geom::CYLINDRICAL: // Laplace operator in cyl coords.
                D_WEIGHT(bi, li) += FDWEIGHT(1, li)*1/xc[around];
                A_WEIGHT(bi, li) += FDWEIGHT(0, li)*1/xc[around];
                break;
            case Geom::SPHERICAL: // Laplace operator in sph coords.
                D_WEIGHT(bi, li) += FDWEIGHT(1, li)*2/xc[around];
                A_WEIGHT(bi, li) += FDWEIGHT(0, li)*2/xc[around];
                break;
            default:
                break;
            }
        }
    }
    delete []c;
}
#undef FDWEIGHT
// D_WEIGHT still defined

#define FACTOR(ri, bi) ( ((ri) < n_factor_affected_k) ? \
            bin_k_factor[bi][i_bin_k[ri]] : 1 )
void
ReactionDiffusion::_fill_local_r(int bi, const double * const __restrict__ y,
                 double * const __restrict__ local_r) const
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
        if (logy)
            local_r[rxni] = exp(local_r[rxni]);

        local_r[rxni] *= FACTOR(rxni, bi)*k[rxni];
    }
}
#undef FACTOR
// D_WEIGHT still defined

// The indices of x, fluxes and bins
// <indices.png>


#define Y(bi, si) y[(bi)*n+(si)]
#define LINC(bi, si) linC[(bi)*n+(si)]
const double *
ReactionDiffusion::_alloc_and_populate_linC(const double * const __restrict__ y) const
{
    int nlinC = n*N;
    double * const linC = (double * const)malloc(nlinC*sizeof(double));
    // TODO: Tune 42...
    ${"#pragma omp parallel for if (N > 42)" if USE_OPENMP else ""}
    for (uint bi=0; bi<N; ++bi)
        for (uint si=0; si<n; ++si)
            LINC(bi, si) = exp(Y(bi, si));
    return linC;
}
// Y, LINC, D_WEIGHT(bi, li) still defined

#define DYDT(bi, si) dydt[(bi)*(n)+(si)]
void
ReactionDiffusion::f(double t, const double * const y, double * const __restrict__ dydt)
{
    // note condifiontal call to free at end of this function
    const double * const linC = (logy && (N > 1)) ? _alloc_and_populate_linC(y) : y;
    if (auto_efield)
        calc_efield(linC);

    ${"double * const local_r = new double[nr];" if not USE_OPENMP else ""}
    ${"#pragma omp parallel for if (N > 2)" if USE_OPENMP else ""}
    for (uint bi=0; bi<N; ++bi){
        // compartment bi
        ${"double * const local_r = new double[nr];" if USE_OPENMP else ""}

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
            // Contributions from diffusion and advection
            // ------------------------------------------
            int starti;
            if ((bi < nsidep) && (!lrefl)){
                starti = 0;
            } else if ((bi >= N-nsidep) && (!rrefl)){
                starti = N - nstencil;
            } else{
                starti = bi - nsidep;
            }
            for (uint si=0; si<n; ++si){ // species index si
                if ((D[si] == 0.0) && ((mobility[si] == 0.0) || efield[bi] == 0.0)) continue;
                double unscaled_diffusion = 0;
                double unscaled_advection = 0;
                for (uint xi=0; xi<nstencil; ++xi){
                    int biw = starti + xi;
                    // reflective logic:
                    if (starti < 0){
                        biw = (biw < 0) ? (-1 - biw) : biw; // lrefl==true
                    } else if (starti >= (int)N - (int)nstencil + 1){
                        biw = (biw >= (int)N) ? (2*N - biw - 1) : biw; // rrefl==true
                    }
                    unscaled_diffusion += D_WEIGHT(bi, xi) * ((logy) ? LINC(biw, si) : Y(biw, si));
                    unscaled_advection += A_WEIGHT(bi, xi) * ((logy) ? LINC(biw, si) : Y(biw, si));
                }
                DYDT(bi, si) += unscaled_diffusion*D[si];
                DYDT(bi, si) += unscaled_advection*efield[bi]*mobility[si];
            }
        }
        if (logy){
            if (logt)
                for (uint si=0; si<n; ++si)
                    DYDT(bi, si) *= exp(t-Y(bi, si));
            else
                for (uint si=0; si<n; ++si)
                    //DYDT(bi, si) *= exp(-Y(bi, si));
                    DYDT(bi, si) /= LINC(bi, si);
        } else {
            if (logt)
                for (uint si=0; si<n; ++si)
                    DYDT(bi, si) *= exp(t);
        }

        ${"delete []local_r;" if USE_OPENMP else ""}

    }
    ${"delete []local_r;" if not USE_OPENMP else ""}

    if (logy && (N > 1))
        free((void*)linC);
}
#undef DYDT
// D_WEIGHT(bi, li), Y(bi, si) and LINC(bi, si) still defined.

#define FOUT(bi, si) fout[(bi)*n+si]
%for token, imaj, imin in [\
    ('dense_jac_rmaj',         '(bri)*n+ri', '(bci)*n + ci'),\
    ('dense_jac_cmaj',         '(bci)*n+ci', '(bri)*n + ri'),\
    ('banded_packed_jac_cmaj', '(bci)*n+ci', '(1+bri-(bci))*n+ri-ci'),\
    ('banded_padded_jac_cmaj', '(bci)*n+ci', '(2+bri-(bci))*n+ri-ci'),\
    ]:
#define JAC(bri, bci, ri, ci) ja[(${imaj})*ldj+${imin}]
void
ReactionDiffusion::${token}(double t, const double * const y,
                            double * const __restrict__ ja, int ldj)
{
    // Note: does not return a strictly correct Jacobian, only 1 pair of bands.
    // `t`: time (log(t) if logt=1)
    // `y`: concentrations (log(conc) if logy=True)
    // `ja`: jacobian (allocated 1D array to hold dense or banded)
    // `ldj`: leading dimension of ja (useful for padding)
    double * const __restrict__ fout = (logy) ? new double[n*N] : nullptr;
    if (fout != nullptr) f(t, y, fout);

    // note condifiontal call to free at end of this function
    const double * const linC = (logy && (N > 1)) ? _alloc_and_populate_linC(y) : y;
    if (auto_efield)
        calc_efield(linC);

    ${'double * const local_r = new double[nr];' if not USE_OPENMP else ''}
    ${'#pragma omp parallel for' if USE_OPENMP else ''}
    for (uint bi=0; bi<N; ++bi){
        // Conc. in `bi:th` compartment
        ${'double * const local_r = new double[nr];' if USE_OPENMP else ''}

        // Contributions from reactions
        // ----------------------------
        _fill_local_r(bi, y, local_r);
        for (uint si=0; si<n; ++si){
            // species si
            for (uint dsi=0; dsi<n; ++dsi){
                // derivative wrt species dsi
                // j_i[si, dsi] = Sum_l(n_lj*Derivative(r[l], local_y[dsi]))
                JAC(bi, bi, si, dsi) = 0.0;
                for (uint rxni=0; rxni<nr; ++rxni){
                    // reaction rxni
                    if (coeff_totl[rxni*n + si] == 0)
                        continue; // species si unaffected by reaction
                    if (coeff_actv[rxni*n + dsi] == 0)
                        continue; // rate of reaction unaffected by species dsi
                    double tmp = coeff_totl[rxni*n + si]*\
                    coeff_actv[rxni*n + dsi]*local_r[rxni];
                    if (!logy)
                        tmp /= Y(bi, dsi);
                    JAC(bi, bi, si, dsi) += tmp;
                }
                if (logy)
                    //JAC(bi, bi, si, dsi) *= exp(-Y(bi, si));
                    JAC(bi, bi, si, dsi) /= LINC(bi, si);
            }
        }

        // Contributions from diffusion
        // ----------------------------
        if (N > 1) {
            uint lbound = _stencil_bi_lbound(bi);
            for (uint si=0; si<n; ++si){ // species index si
                if ((D[si] == 0.0) && ((mobility[si] == 0.0) || efield[bi] == 0.0)) continue;
                // Not a strict Jacobian only block diagonal plus closest bands...
                for (uint k=0; k<nstencil; ++k){
                    const uint sbi = _xc_bi_map(lbound+k);
                    if (sbi == bi) {
                        JAC(bi, bi, si, si) += D[si]*D_WEIGHT(bi, k);
                        JAC(bi, bi, si, si) += efield[bi]*mobility[si]*A_WEIGHT(bi, k);
                    } else {
                        if (bi > 0)
                            if (sbi == bi-1){
                                double Cfactor = (logy ? LINC(bi-1, si)/LINC(bi, si) : 1.0);
                                JAC(bi, bi-1, si, si) += D[si]*D_WEIGHT(bi, k)*Cfactor;
                                JAC(bi, bi-1, si, si) += efield[bi]*mobility[si]*A_WEIGHT(bi, k)*\
                                    Cfactor;
                            }
                        if (bi < N-1)
                            if (sbi == bi+1){
                                double Cfactor = (logy ? LINC(bi+1, si)/LINC(bi, si) : 1.0);
                                JAC(bi, bi+1, si, si) += D[si]*D_WEIGHT(bi, k)*Cfactor;
                                JAC(bi, bi+1, si, si) += efield[bi]*mobility[si]*A_WEIGHT(bi, k)*\
                                    Cfactor; 
                            }
                    }
                }
            }
        }

        // Logartihmic time
        // ----------------------------
        if (logt || logy){
            for (uint si=0; si<n; ++si){
                if (logt){
                    for (uint dsi=0; dsi<n; ++dsi)
                        JAC(bi, bi, si, dsi) *= exp(t);
                    if (bi>0)
                        JAC(bi, bi - 1, si, si) *= exp(t);
                    if (bi<N-1)
                        JAC(bi, bi + 1, si, si) *= exp(t);
                }
                if (logy)
                    JAC(bi, bi, si, si) -= FOUT(bi, si);
            }
        }
        ${'delete []local_r;' if USE_OPENMP else ''}
    }
    ${'delete []local_r;' if not USE_OPENMP else ''}
    if (logy)
        delete []fout;
    if (logy && (N > 1))
        free((void*)linC);
}
#undef JAC
%endfor

#undef FOUT
#undef LINC
#undef Y
#undef D_WEIGHT
#undef A_WEIGHT

void ReactionDiffusion::per_rxn_contrib_to_fi(double t, const double * const __restrict__ y,
                                              uint si, double * const __restrict__ out) const
{
    double * const local_r = new double[nr];
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

void ReactionDiffusion::calc_efield(const double * const linC)
{
    // Prototype for self-generated electric field
    double netchg[N];
    const double F = 96485.3399; // Faraday's constant, [C/mol]
    const double pi = 3.14159265358979324;
    double Q = surf_chg;
    double nx, cx = logx ? exp(x[0]) : x[0];
    for (uint bi=0; bi<N; ++bi){
        netchg[bi] = 0.0;
        for (uint si=0; si<n; ++si)
            netchg[bi] += z_chg[si]*linC[bi*n+si];
    }
    for (uint bi=0; bi<N; ++bi){
        const double r = logx ? exp(xc[nsidep+bi]) : xc[nsidep+bi];
        nx = logx ? exp(x[bi+1]) : x[bi+1];
        switch(geom){
        case Geom::FLAT:
            efield[bi] = F*Q;
            Q += netchg[bi]*(nx - cx);
            break;
        case Geom::CYLINDRICAL:
            efield[bi] = F*Q/(2*pi*eps*r); // Gauss's law
            Q += netchg[bi]*pi*(nx*nx - cx*cx);
            break;
        case Geom::SPHERICAL:
            efield[bi] = F*Q/(4*pi*eps*r*r); // Gauss's law
            Q += netchg[bi]*4*pi/3*(nx*nx*nx - cx*cx*cx);
            break;
        }
        cx = nx;
    }
    if (geom == Geom::FLAT){
        Q = 0.0;
        for (uint bi=N; bi>0; --bi){ // unsigned int..
            nx = logx ? exp(x[bi-1]) : x[bi-1];
            efield[bi-1] += F*Q;
            Q += netchg[bi]*(cx - nx);
            cx = nx;
        }
    }
}

}; // namespace chemreac
