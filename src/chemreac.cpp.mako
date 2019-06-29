// ${'-{0}- eval: ({1}) -{0}-'.format('*', 'read-only-mode')}
// ${__import__('codecs').encode('Guvf svyr jnf trarengrq, qb abg rqvg', 'rot_13')}
<%doc> This is a source file template for use with the Python rendering engine "mako" </%doc>
#include <algorithm> // std::count
//#include <vector>    // std::vector
#include <algorithm> // std::max, std::min
#include <cstdlib> // free,  C++11 aligned_alloc
#include <memory>
#include "anyode/anyode_buffer.hpp"
#include "anyode/anyode_decomposition.hpp"
#include "finitediff_templated.hpp" // fintie differences
#include "chemreac.hpp"

#include <iostream> //DEBUG


#if defined(CHEMREAC_WITH_DATA_DUMPING)
#include <cstdio>
#include <iostream>
#include <sstream>
#include <iomanip>
#define PRINT_ARR(ARR, LARR) for(int i_=0; i_<LARR; ++i_) {std::cout << ARR[i_] << " ";}; std::cout << '\n';
#include "chemreac_util.hpp" // save_array, load_array
#endif

%if WITH_OPENMP:
#include <omp.h>
%else:
#define omp_get_thread_num() 0
%endif


namespace chemreac {
using std::vector;
using std::count;
using std::min;
using std::max;

#include <cstdio>

#define expb(arg) (use_log2 ? std::exp2(arg) : std::exp(arg))
#define logb(arg) (use_log2 ? std::log2(arg) : std::pow(2, arg))

#define GRAD_WEIGHT(bi, li) grad_weight[nstencil*(bi) + li]


// 1D discretized reaction diffusion
template<typename Real_t>
ReactionDiffusion<Real_t>::ReactionDiffusion(
    int n,
    const vector<vector<int> > stoich_active,
    const vector<vector<int> > stoich_prod,
    vector<Real_t> k,
    int N,
    vector<Real_t> D,
    const vector<int> z_chg,
    vector<Real_t> mobility,
    const vector<Real_t> x, // separation
    vector<vector<int> > stoich_inact, // vectors of size 0 in stoich_actv_ => "copy from stoich_reac"
    int geom_,
    bool logy,
    bool logt,
    bool logx,
    int nstencil,
    bool lrefl,
    bool rrefl,
    bool auto_efield,
    pair<Real_t, Real_t> surf_chg,
    Real_t eps_rel,
    Real_t faraday_const,
    Real_t vacuum_permittivity,
    vector<vector<Real_t>> g_values,
    vector<int> g_value_parents,
    vector<vector<Real_t>> fields,
    vector<int> modulated_rxns,
    vector<vector<Real_t> > modulation,
    Real_t ilu_limit,
    int n_jac_diags,
    bool use_log2,
    bool clip_to_pos):
    n(n), N(N), nstencil(nstencil), nsidep((nstencil-1)/2), nr(stoich_active.size()),
    coeff_active(buffer_factory<int>(nr*n)),
    coeff_prod(buffer_factory<int>(nr*n)),
    coeff_total(buffer_factory<int>(nr*n)),
    coeff_inact(buffer_factory<int>(nr*n)),
    lap_weight(buffer_factory<double>(nstencil*N)),
    div_weight(buffer_factory<double>(nstencil*N)),
    grad_weight(buffer_factory<double>(nstencil*N)),
    efield(buffer_factory<double>(N)),
    netchg(buffer_factory<double>(N)),
    xc(buffer_factory<double>(nsidep + N + nsidep)),
    work1(buffer_factory<double>(n*N)),
    work2(buffer_factory<double>(n*N)),
    work3(buffer_factory<double>(${"((nr/8)+1)*8*omp_get_num_threads()" if WITH_OPENMP else "nr"})),
    logy(logy), logt(logt), logx(logx), stoich_active(stoich_active),
    stoich_inact(stoich_inact), stoich_prod(stoich_prod),
    k(k),  D(D), z_chg(z_chg), mobility(mobility), x(x), lrefl(lrefl), rrefl(rrefl),
    auto_efield(auto_efield),
    surf_chg(surf_chg), eps_rel(eps_rel), faraday_const(faraday_const),
    vacuum_permittivity(vacuum_permittivity),
    g_value_parents(g_value_parents), modulated_rxns(modulated_rxns), modulation(modulation),
    ilu_limit(ilu_limit), n_jac_diags((n_jac_diags == 0) ? nsidep : n_jac_diags), use_log2(use_log2),
    clip_to_pos(clip_to_pos)
{
    if (N < 1) throw std::logic_error("One is the smallest number of bins.");
    if (N == 2) throw std::logic_error("2nd order PDE requires at least 3 stencil points.");
    if (nstencil % 2 == 0) throw std::logic_error("Only odd number of stencil points supported");
    if ((N == 1) && (nstencil != 1)) throw std::logic_error("You must set nstencil=1 for N=1");
    if ((N > 1) && (nstencil <= 1)) throw std::logic_error("You must set nstencil>1 for N>1");
    if (stoich_active.size() != stoich_prod.size())
        throw std::length_error(
            "stoich_active and stoich_prod of different sizes.");
    if (k.size() != stoich_prod.size())
        throw std::length_error(
            "k and stoich_prod of different sizes.");
    if (N>1){
        if (D.size() != (unsigned)(n*N))
            throw std::length_error(
                "Length of D does not match number of species * number of bins.");
        if (mobility.size() != (unsigned)n)
            throw std::length_error(
                "Length of mobility does not match number of species.");
        if (z_chg.size() != (unsigned)n)
            throw std::length_error(
                "Length of z_chg does not match number of species.");
        if (x.size() != (unsigned)(N + 1))
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
    case 3:
        geom = Geom::PERIODIC;
        break;
    default:
        throw std::logic_error("Unknown geom.");
    }

    // Finite difference scheme
    for (int i=0; i<N; ++i) efield[i] = 0.0;
    for (int i=0; i<N; ++i)
        xc[nsidep + i] = (x[i] + x[i + 1])/2;

    for (int i=0; i<nsidep; ++i){
        // reflection
        xc[nsidep - i - 1] = 2*x[0] - xc[nsidep + i];
        xc[nsidep + i + N] = 2*x[N] - xc[nsidep + N - i - 1];
    }

    // Precalc coeffs for Jacobian for current geom.
    // not centered diffs close to boundaries
    for (int bi=0; bi<N; bi++)
        apply_fd_(bi);

    gradD = buffer_factory<Real_t>(N*n);
    for (int bi=0; bi<N; ++bi){
        int starti = start_idx_(bi);
        for (int si=0; si<n; ++si){
            gradD[bi*n + si] = 0;
            for (int li=0; li<nstencil; ++li){
                int biw = biw_(starti, li);
                gradD[bi*n + si] += GRAD_WEIGHT(bi, li)*D[biw*n + si];
            }
        }
    }

    // Stoichiometry
    for (int ri=0; ri<nr; ++ri){
        for (auto si=stoich_active[ri].begin(); si != stoich_active[ri].end(); ++si)
            if (*si > n-1)
                throw std::logic_error("At least one species index in stoich_active > (n-1)");
        for (auto si=stoich_prod[ri].begin(); si != stoich_prod[ri].end(); ++si)
            if (*si > n-1)
                throw std::logic_error("At least one species index in stoich_prod > (n-1)");
        for (auto si=stoich_inact[ri].begin(); si != stoich_inact[ri].end(); ++si)
            if (*si > n-1)
                throw std::logic_error("At least one species index in stoich_inact > (n-1)");
    }

    stoich_inact.reserve(nr);
    for (int rxni=0; rxni<nr; ++rxni){ // reaction index
        for (int si=0; si<n; ++si){ // species index
            coeff_active[rxni*n+si] = count(stoich_active[rxni].begin(), stoich_active[rxni].end(), si);
            coeff_inact[rxni*n+si] = count(stoich_inact[rxni].begin(), stoich_inact[rxni].end(), si);
            coeff_prod[rxni*n+si] = count(stoich_prod[rxni].begin(), stoich_prod[rxni].end(), si);
            coeff_total[rxni*n+si] = coeff_prod[rxni*n+si] - coeff_active[rxni*n+si] - coeff_inact[rxni*n+si];
        }
    }

    // Handle g_values
    if (fields.size() != g_values.size())
        throw std::logic_error("fields and g_values need to be of equal length");
    if (this->g_value_parents.size() != g_values.size())
        throw std::logic_error("g_value_parents and g_values need to be of equal length");

    for (const auto& gs : g_values)
        if (gs.size() != (unsigned)n)
            throw std::logic_error("vectors in g_values need to be of length n");

    for (const auto& fs : fields)
        if (fs.size() != (unsigned)N)
            throw std::logic_error("vectors in fields need to be of length N");

    this->g_values = g_values;
    this->fields = fields;

    // Sanity check modulation
    for (const auto rxni : this->modulated_rxns)
        if (rxni >= (int)nr || rxni < 0)
            throw std::logic_error("illegal reaction index in modulated_rxns");
    if (this->modulation.size() != this->modulated_rxns.size())
        throw std::logic_error("modulation size differs from modulated_rxns");
    for (const auto& mdltn : this->modulation)
        if (mdltn.size() != (unsigned)N)
            throw std::logic_error("illegally sized vector in modulation");
}

template<typename Real_t>
ReactionDiffusion<Real_t>::~ReactionDiffusion(){}

template<typename Real_t>
int ReactionDiffusion<Real_t>::start_idx_(int bi) const {
    int starti;
    if ((bi < this->nsidep) && (!this->lrefl)){
        starti = 0;
    } else if ((bi >= this->N - this->nsidep) && (!this->rrefl)){
        starti = this->N - this->nstencil;
    } else{
        starti = bi - this->nsidep;
    }
    return starti;
}

template<typename Real_t>
int ReactionDiffusion<Real_t>::biw_(int starti, int li) const {
    int biw = starti + li;
    // reflective logic:
    if (starti < 0){
        biw = (biw < 0) ? (-1 - biw) : biw; // lrefl==true
    } else if (starti >= (int)N - (int)nstencil + 1){
        biw = (biw >= (int)N) ? (2*N - biw - 1) : biw; // rrefl==true
    }
    return biw;
}

template<typename Real_t>
void
ReactionDiffusion<Real_t>::zero_counters(){
    nfev = 0;
    njev = 0;
    nprec_setup = 0;
    nprec_solve = 0;
    njacvec_dot = 0;
    nprec_solve_ilu = 0;
    nprec_solve_lu = 0;
}

template<typename Real_t>
int
ReactionDiffusion<Real_t>::get_ny() const
{
    return n*N;
}

template<typename Real_t>
int
ReactionDiffusion<Real_t>::get_mlower() const
{
    if (N > 1)
        return n*n_jac_diags;
    else
        return -1;
}

template<typename Real_t>
int
ReactionDiffusion<Real_t>::get_mupper() const
{
    return this->get_mlower();
}

template<typename Real_t>
Real_t
ReactionDiffusion<Real_t>::get_dx_max(Real_t x, const Real_t * const y)
{
    %for side in "upper lower".split():
    if (m_${side}_bounds.size() != static_cast<std::size_t>(n*N)) {
        throw std::runtime_error("trying to use get_dx_max without m_${side}_bounds set to correct size.");
    }
    %endfor
    const int ny = get_ny();
    auto fvec = std::vector<Real_t>(ny);
    auto hvec = std::vector<Real_t>(ny);
    rhs(x, y, &fvec[0]);
    for (int idx=0; idx < ny; ++idx){
        if (fvec[idx] == 0) {
            hvec[idx] = std::numeric_limits<Real_t>::max();
        } else if (fvec[idx] > 0) {
            hvec[idx] = std::abs(m_upper_bounds[idx] - y[idx])/fvec[idx];
        } else { // fvec[idx] < 0
            hvec[idx] = std::abs((m_lower_bounds[idx] - y[idx])/fvec[idx]);
        }
    }
    const auto result = *std::min_element(std::begin(hvec), std::end(hvec));
    if (m_get_dx_max_factor == 0.0)
        return result;
    else if (m_get_dx_max_factor < 0.0)
        return -m_get_dx_max_factor*result;
    else
        return m_get_dx_max_factor*result;
}

template<typename Real_t>
Real_t
ReactionDiffusion<Real_t>::get_dx0(Real_t x, const Real_t * const y)
{
    %for side in "upper lower".split():
    if (m_${side}_bounds.size() != static_cast<std::size_t>(n*N)) {
        return this->default_dx0;
    }
    %endfor
    return m_get_dx0_factor*std::min(get_dx_max(x, y), m_get_dx0_max_dx);
}

template<typename Real_t>
int
ReactionDiffusion<Real_t>::stencil_bi_lbound_(int bi) const
{
    const int le = lrefl ? 0 : nsidep;
    const int re = rrefl ? 0 : nsidep;
    return max(le, min(N + 2*nsidep - re - nstencil, bi));
}

template<typename Real_t>
int
ReactionDiffusion<Real_t>::xc_bi_map_(int xci) const
{
    if (xci < nsidep)
        return nsidep - xci - 1;
    else if (xci >= N+nsidep)
        return nsidep + 2*N - xci - 1;
    else
        return xci - nsidep;
}


#define LAP_WEIGHT(bi, li) lap_weight[nstencil*(bi) + li]
#define DIV_WEIGHT(bi, li) div_weight[nstencil*(bi) + li]
#define FDWEIGHT(order, local_index) c[nstencil*(order) + local_index]
template<typename Real_t>
void
ReactionDiffusion<Real_t>::apply_fd_(int bi){
    auto c = AnyODE::buffer_factory<Real_t>(3*nstencil);
    auto lxc = AnyODE::buffer_factory<Real_t>(nstencil); // local shifted x-centers
    int around = bi + nsidep;
    int start = bi;
    if (!lrefl) // shifted finite diff
        start = max(nsidep, start);
    if (!rrefl) // shifted finite diff
        start = min(N - nstencil + nsidep, start);

    for (int li=0; li<nstencil; ++li) // li: local index
        lxc[li] = xc[start + li] - xc[around];
    finitediff::populate_weights<Real_t>(0, AnyODE::buffer_get_raw_ptr(lxc), nstencil-1, 2, AnyODE::buffer_get_raw_ptr(c));

    const Real_t logbdenom = use_log2 ? 1/log(2) : 1;

    for (int li=0; li<nstencil; ++li){ // li: local index
        LAP_WEIGHT(bi, li) = FDWEIGHT(2, li);
        DIV_WEIGHT(bi, li) = FDWEIGHT(1, li);
        GRAD_WEIGHT(bi, li) = FDWEIGHT(1, li); // == DIV_WEIGHT if geom == 'f'
        if (logx){
            LAP_WEIGHT(bi, li) *= logbdenom*logbdenom;
            DIV_WEIGHT(bi, li) *= logbdenom;
            GRAD_WEIGHT(bi, li) *= logbdenom;
            switch(geom){
            case Geom::FLAT:
            case Geom::PERIODIC:
                LAP_WEIGHT(bi, li) -= FDWEIGHT(1, li)*logbdenom;
                break;
            case Geom::CYLINDRICAL:
                DIV_WEIGHT(bi, li) += FDWEIGHT(0, li);
                break;
            case Geom::SPHERICAL:
                LAP_WEIGHT(bi, li) += FDWEIGHT(1, li)*logbdenom;
                DIV_WEIGHT(bi, li) += 2*FDWEIGHT(0, li);
                break;
            }
            LAP_WEIGHT(bi, li) *= expb(-2*xc[around]);
            DIV_WEIGHT(bi, li) *= expb(-xc[around]);
            GRAD_WEIGHT(bi, li) *= expb(-xc[around]);
        } else {
            switch(geom){
            case Geom::CYLINDRICAL: // Laplace operator in cyl coords.
                LAP_WEIGHT(bi, li) += FDWEIGHT(1, li)*1/xc[around];
                DIV_WEIGHT(bi, li) += FDWEIGHT(0, li)*1/xc[around];
                break;
            case Geom::SPHERICAL: // Laplace operator in sph coords.
                LAP_WEIGHT(bi, li) += FDWEIGHT(1, li)*2/xc[around];
                DIV_WEIGHT(bi, li) += FDWEIGHT(0, li)*2/xc[around];
                break;
            default:
                break;
            }
        }
    }
}
#undef FDWEIGHT

template<typename Real_t>
Real_t
ReactionDiffusion<Real_t>::get_mod_k(int bi, int ri) const{
    Real_t tmp = k[ri];
    // Modulation
    int enumer = -1;
    for (auto mi : this->modulated_rxns){
        enumer++;
        if (mi == (int)ri)
            tmp *= this->modulation[enumer][bi];
    }
    return tmp;
}

template<typename Real_t>
void
ReactionDiffusion<Real_t>::fill_local_r_(int bi, const Real_t * const ANYODE_RESTRICT C,
                                         Real_t * const ANYODE_RESTRICT local_r) const
{
    // intent(out) :: local_r
    for (int rxni=0; rxni<nr; ++rxni){
        // reaction rxni
        Real_t tmp = 1;

        // Kinetically active reactants (law of massaction)
        for (unsigned rnti=0; rnti < stoich_active[rxni].size(); ++rnti){
            // reactant index rnti
            int si = stoich_active[rxni][rnti];
            tmp *= C[bi*n+si];
        }
        // Rate constant
        local_r[rxni] = get_mod_k(bi, rxni)*tmp;
    }
}

// The indices of x, fluxes and bins
// <indices.png>


#define Y(bi, si) y[(bi)*n+(si)]
template<typename Real_t>
void
ReactionDiffusion<Real_t>::populate_linC(Real_t * const ANYODE_RESTRICT linC,
                                         const Real_t * const ANYODE_RESTRICT y,
                                         bool apply_exp, bool recip) const
{
    ${"#pragma omp parallel for schedule(static) if (N*n > 65536)" if WITH_OPENMP else ""}
    for (int bi=0; bi<N; ++bi){
        if (recip) {
            if (apply_exp) {
                for (int si=0; si<n; ++si){
                    linC[bi*n + si] = expb(-Y(bi, si));
                }
            } else {
                for (int si=0; si<n; ++si){
                    linC[bi*n + si] = 1.0/Y(bi, si);
                }
            }
        } else {
            if (apply_exp) {
                for (int si=0; si<n; ++si){
                    linC[bi*n + si] = expb(Y(bi, si));
                }
            } else {
                for (int si=0; si<n; ++si){
                    linC[bi*n + si] = Y(bi, si);
                }
            }
        }
    }
}
#define LINC(bi, si) linC[(bi)*n+(si)]
#define RLINC(bi, si) rlinC[(bi)*n+(si)]

#define DYDT(bi, si) dydt[(bi)*(n)+(si)]
template<typename Real_t>
AnyODE::Status
ReactionDiffusion<Real_t>::rhs(Real_t t, const Real_t * const y, Real_t * const ANYODE_RESTRICT dydt)
{
    // note condifiontal call to free at end of this function
    bool use_work = false;
    if (logy) {
        populate_linC(AnyODE::buffer_get_raw_ptr(work1), y, true);
        populate_linC(AnyODE::buffer_get_raw_ptr(work2), y, true, true);
        use_work = true;
    }
    if (clip_to_pos) {
        if (!use_work) {
            memcpy(AnyODE::buffer_get_raw_ptr(work1), y, sizeof(Real_t)*n*N);
            for (int i=0; i<get_ny(); ++i){
                work2[i] = 1.0/y[i];
            }
            use_work = true;
        }
        for (int i=0; i<get_ny(); ++i){
            work1[i] = (y[i] < 0) ? 0 : y[i];
        }
        for (int i=0; i<get_ny(); ++i){
            work2[i] = (y[i] < 0) ? INFINITY : y[i];
        }
    }
    const Real_t * const linC = (use_work) ? AnyODE::buffer_get_raw_ptr(work1) : y;
    const Real_t * const rlinC = (use_work) ? AnyODE::buffer_get_raw_ptr(work2) : nullptr;
    if (m_error_outside_bounds) {
        if (m_lower_bounds.size() > 0) {
            for (int i=0; i < n*N; ++i) {
                if (linC[i] < m_lower_bounds[i]) {
                    std::clog << "Lower bound (" << m_lower_bounds[0] << ") for "
                              << std::to_string(i)
                              << " not satisfied (" << linC[i] << ") at t="<< t << "\n";
                    return AnyODE::Status::recoverable_error;
                }
            }
        }
        if (m_upper_bounds.size() > 0) {
            for (int i=0; i < n*N; ++i) {
                if (linC[i] > m_upper_bounds[i]) {
                    std::clog << "Upper bound (" << m_upper_bounds[0] << ") for "
                              <<  std::to_string(i)
                              << " not satisfied (" << linC[i] << ") at t="<< t << "\n";
                    return AnyODE::Status::recoverable_error;
                }
            }
        }
    }
    if (auto_efield){
        calc_efield(linC);
    }
    const Real_t expb_t = (logt) ? expb(t) : 0.0;
    ${"Real_t * const local_r = AnyODE::buffer_get_raw_ptr(work3);" if not WITH_OPENMP else ""}
    ${"#pragma omp parallel for schedule(static) if (N*n > 65536)" if WITH_OPENMP else ""}
    for (int bi=0; bi<N; ++bi){
        // compartment bi
        ${"Real_t * const local_r = AnyODE::buffer_get_raw_ptr(work3) + ((nr/8)+1)*8*omp_get_thread_num();" if WITH_OPENMP else ""}

        for (int si=0; si<n; ++si)
            DYDT(bi, si) = 0.0; // zero out

        // Contributions from reactions
        // ----------------------------
        fill_local_r_(bi, linC, local_r);
        for (int rxni=0; rxni<nr; ++rxni){
            // reaction index rxni
            for (int si=0; si<n; ++si){
                // species index si
                const int overall = coeff_total[rxni*n + si];
                if (overall != 0)
                    DYDT(bi, si) += overall*local_r[rxni];
            }
        }
        // Contribution from particle/electromagnetic fields
        for (unsigned fi=0; fi<this->fields.size(); ++fi){
            if (fields[fi][bi] == 0)
                continue; // exit early
            const Real_t gfact = (g_value_parents[fi] == -1) ? \
                1.0 : LINC(bi, g_value_parents[fi]);
            for (int si=0; si<n; ++si)
                if (g_values[fi][si] != 0)
                    DYDT(bi, si) += fields[fi][bi]*g_values[fi][si]*gfact;
        }

        if (N>1){
            // Contributions from diffusion and advection
            // ------------------------------------------
            int starti = start_idx_(bi);
            for (int si=0; si<n; ++si){ // species index si
                if ((D[bi*n + si] == 0.0) && (mobility[si] == 0.0)) continue;
                Real_t diffusion_unscaled = 0;
                Real_t diffusion_correction = 0;
                Real_t advection_unscaled = 0;
                for (int xi=0; xi<nstencil; ++xi){
                    int biw = biw_(starti, xi);
                    diffusion_unscaled += LAP_WEIGHT(bi, xi) * LINC(biw, si);
                    diffusion_correction += GRAD_WEIGHT(bi, xi) * LINC(biw, si);
                    advection_unscaled += DIV_WEIGHT(bi, xi) * \
                        (LINC(biw, si)*efield[bi] + LINC(bi, si)*efield[biw]);
                }
                DYDT(bi, si) += diffusion_unscaled*D[bi*n + si];
                DYDT(bi, si) += diffusion_correction*gradD[bi*n + si];
                DYDT(bi, si) += advection_unscaled*-mobility[si];
            }
        }
        for (int si=0; si<n; ++si){
            if (logy){
                DYDT(bi, si) *= RLINC(bi, si);
                if (!logt and use_log2)
                    DYDT(bi, si) /= log(2);
            }
            if (logt){
                DYDT(bi, si) *= expb_t;
                if (!logy and use_log2)
                    DYDT(bi, si) *= log(2);
            }
        }
    }
    nfev++;
    return AnyODE::Status::success;
}
#undef DYDT

#define FOUT(bi, si) fout[(bi)*n+si]
#define SUP(di, bi, li) jac.sup(di, bi, li)
%for token in ["dense_jac_rmaj", "dense_jac_cmaj", "banded_jac_cmaj", "compressed_jac_cmaj"]:
template<typename Real_t>
AnyODE::Status
ReactionDiffusion<Real_t>::${token}(Real_t t,
                                    const Real_t * const ANYODE_RESTRICT y,
                                    const Real_t * const ANYODE_RESTRICT fy,
                                    Real_t * const ANYODE_RESTRICT ja, long int ldj
                                    ${', double * const ANYODE_RESTRICT /* dfdt */' if token.startswith('dense') else ''})
{
    // Note: blocks are zeroed out, diagonals only incremented
    // `t`: time (log(t) if logt=1)
    // `y`: concentrations (log(conc) if logy=True)
    // `ja`: jacobian (allocated 1D array to hold dense or banded)
    // `ldj`: leading dimension of ja (useful for padding, ignored by compressed_*)
 %if token.startswith("compressed"):
    ignore(ldj);
    const int nsat = (geom == Geom::PERIODIC) ? nsidep : 0 ;
    const int ld = n;
    block_diag_ilu::BlockDiagMatrix<Real_t> jac {ja, N, n, n_jac_diags, nsat, ld};
 %elif token.startswith("banded_jac_cmaj"):
    block_diag_ilu::BlockBandedMatrix<Real_t> jac {ja-n*n_jac_diags, N, n, n_jac_diags, static_cast<int>(ldj)};
    // no need to zero out jac (see documentation for CVDlsJacFn)
    //std::memset(jac - n*n_jac_diags, 0, sizeof(Real_t)*ldj*N*n);
    //jac.set_to(0.0);
 %elif token.startswith("dense_jac_cmaj"):
    block_diag_ilu::BlockDenseMatrix<Real_t> jac {ja, N, n, n_jac_diags, static_cast<int>(ldj), true};
 %elif token.startswith("dense_jac_rmaj"):
    block_diag_ilu::BlockDenseMatrix<Real_t> jac {ja, N, n, n_jac_diags, static_cast<int>(ldj), false};
 %else:
    #error "Unhandled token."
 %endif
    const Real_t exp_t = (logt) ? expb(t) : 0.0;
    const Real_t logbfactor = use_log2 ? log(2) : 1;

    Real_t * fout = nullptr;
    if (logy){ // fy useful..
        if (fy){
            fout = const_cast<Real_t *>(fy);
        } else {
            fout = new Real_t[n*N];
            rhs(t, y, fout);
        }
    }
    bool use_work = false;
    // note conditional call to free at end of this function
    if (logy) {
        populate_linC(AnyODE::buffer_get_raw_ptr(work1), y, true, false);
        populate_linC(AnyODE::buffer_get_raw_ptr(work2), y, true, true);
        use_work = true;
    }
    if (clip_to_pos) {
        if (!use_work) {
            memcpy(AnyODE::buffer_get_raw_ptr(work1), y, sizeof(Real_t)*n*N);
            for (int i=0; i<get_ny(); ++i){
                work2[i] = 1.0/y[i];
            }
            use_work = true;
        }
        for (int i=0; i<get_ny(); ++i){
            work1[i] = (y[i] < 0) ? 0 : y[i];
        }
        for (int i=0; i<get_ny(); ++i){
            work2[i] = (y[i] < 0) ? INFINITY : y[i];
        }
    }
    const Real_t * const linC = (use_work) ? AnyODE::buffer_get_raw_ptr(work1) : y;
    const Real_t * const rlinC = (use_work) ? AnyODE::buffer_get_raw_ptr(work2) : nullptr;
    if (auto_efield) {
        calc_efield(linC);
    }

    ${"#pragma omp parallel for schedule(static) if (N*n*n > 65536)" if WITH_OPENMP else ""}
    for (int bi=0; bi<N; ++bi){
        // Conc. in `bi:th` compartment
        // Contributions from reactions and fields
        // ---------------------------------------
        for (int si=0; si<n; ++si){
            // species si
            for (int dsi=0; dsi<n; ++dsi){
                // derivative wrt species dsi
                jac.block(bi, si, dsi) = 0.0;
                for (int rxni=0; rxni<nr; ++rxni){
                    // reaction rxni
                    const int Akj = coeff_active[rxni*n + dsi];
                    const int Ski = coeff_total[rxni*n + si];
                    if (Akj == 0 || Ski == 0)
                        continue;
                    Real_t qkj = get_mod_k(bi, rxni)*Akj*pow(LINC(bi, dsi), Akj-1);
                    for (unsigned rnti=0; rnti < stoich_active[rxni].size(); ++rnti){
                        const int rnti_si = stoich_active[rxni][rnti];
                        if (rnti_si == dsi)
                            continue;
                        qkj *= LINC(bi, rnti_si);
                    }
                    jac.block(bi, si, dsi) += Ski*qkj;
                    //std::cout << "jac.block(" << bi << ", " << si << ", " << dsi <<") = " << jac.block(bi, si, dsi) << "\n";
                }
                // Contribution from particle/electric fields
                for (unsigned fi=0; fi<(this->fields.size()); ++fi){
                    const int Ski = (g_values[fi][si] != 0.0) ? 1 : 0;
                    const int Akj = ((int)dsi == g_value_parents[fi]) ? 1 : 0;
                    const Real_t rk = fields[fi][bi]*g_values[fi][si];
                    if (Ski == 0 || Akj == 0 || rk == 0)
                        continue;
                    jac.block(bi, si, dsi) += Akj*Ski*rk;
                }
            }
        }


        // Contributions from diffusion
        // ----------------------------
        if (N > 1) {
            int lbound = stencil_bi_lbound_(bi);
            for (int si=0; si<n; ++si){ // species index si
                if ((D[bi*n + si] == 0.0) && (mobility[si] == 0.0)) continue; // exit early if possible
                for (int k=0; k<nstencil; ++k){
                    const int sbi = xc_bi_map_(lbound+k);
                    jac.block(bi, si, si) += -mobility[si]*efield[sbi]*DIV_WEIGHT(bi, k);
                    if (sbi == bi) {
                        jac.block(bi, si, si) += D[bi*n + si]*LAP_WEIGHT(bi, k);
                        jac.block(bi, si, si) += gradD[bi*n + si]*GRAD_WEIGHT(bi, k);
                        jac.block(bi, si, si) += efield[bi]*-mobility[si]*DIV_WEIGHT(bi, k);
                    } else {
                        for (int di=0; di<n_jac_diags; ++di){
                            if ((bi >= di+1) and (sbi == bi-di-1)){
                                jac.sub(di, bi-di-1, si) += D[bi*n + si]*LAP_WEIGHT(bi, k);
                                jac.sub(di, bi-di-1, si) += gradD[bi*n+si]*GRAD_WEIGHT(bi, k);
                                jac.sub(di, bi-di-1, si) += efield[bi]*-mobility[si]*DIV_WEIGHT(bi, k);
                            }
                            if ((bi < N-di-1) and (sbi == bi+di+1)){
                                jac.sup(di, bi, si) += D[bi*n + si]*LAP_WEIGHT(bi, k);
                                jac.sup(di, bi, si) += gradD[bi*n + si]*GRAD_WEIGHT(bi, k);
                                jac.sup(di, bi, si) += efield[bi]*-mobility[si]*DIV_WEIGHT(bi, k);
                            }
                        }
                    }
                }
            }
        }

        // Logartihmic transformations
        // ---------------------------
        if (logy || logt){
            for (int si=0; si<n; ++si){
                for (int dsi=0; dsi<n; ++dsi){
                    if (logy){
                        jac.block(bi, si, dsi) *= LINC(bi, dsi)*RLINC(bi, si);
                    }
                    if (logt)
                        jac.block(bi, si, dsi) *= exp_t*logbfactor;
                    if (logy && dsi == si)
                        jac.block(bi, si, si) -= FOUT(bi, si)*logbfactor;
                }
                for (int di=0; di<n_jac_diags; ++di){
                    if (bi > di){
                        if (logy)
                            jac.sub(di, bi-di-1, si) *= LINC(bi-di-1, si)*RLINC(bi, si);
                        if (logt)
                            jac.sub(di, bi-di-1, si) *= exp_t*logbfactor;
                    }
                    if (bi < N-di-1){
                        if (logy)
                            jac.sup(di, bi, si) *= LINC(bi+di+1, si)*RLINC(bi, si);
                        if (logt)
                            jac.sup(di, bi, si) *= exp_t*logbfactor;
                    }
                }
            }
        }
    }
    if (logy && !fy)
        delete []fout;
    njev++;
#if defined(CHEMREAC_WITH_DATA_DUMPING)
    std::ostringstream fname;
    fname << "jac_" << std::setfill('0') << std::setw(5) << njev << ".dat";
  %if token.startswith("banded_jac_cmaj"):
    save_array(jac.m_data, ldj*n*N, fname.str());
  %else:
    save_array(ja, ldj*n*N, fname.str());
  %endif
#endif
    return AnyODE::Status::success;
}
%endfor
#undef FOUT


template<typename Real_t>
AnyODE::Status
ReactionDiffusion<Real_t>::jac_times_vec(const Real_t * const ANYODE_RESTRICT vec,
                                         Real_t * const ANYODE_RESTRICT out,
                                         Real_t t,
                                         const Real_t * const ANYODE_RESTRICT y,
                                         const Real_t * const ANYODE_RESTRICT fy
                                         )
{
    // See 4.6.7 on page 67 (77) in cvs_guide.pdf (Sundials 2.5)
    ignore(t);
    if (!jac_times_cache){
        const int nsat = (geom == Geom::PERIODIC) ? nsidep : 0;
        const int ld = n;
        jac_times_cache = AnyODE::make_unique<block_diag_ilu::BlockDiagMatrix<Real_t>>(nullptr, N, n, nsidep, nsat, ld);
        jac_times_cache->set_to(0.0); // compressed_jac_cmaj only increments diagonals
        const int ld_dummy = 0;
        compressed_jac_cmaj(t, y, fy, jac_times_cache->m_data, ld_dummy);
    }
    std::memset(out, 0, sizeof(Real_t)*get_ny());
    jac_times_cache->dot_vec(vec, out);
    njacvec_dot++;
    return AnyODE::Status::success;
}

template<typename Real_t>
AnyODE::Status
ReactionDiffusion<Real_t>::prec_setup(Real_t t,
                                      const Real_t * const ANYODE_RESTRICT y,
                                      const Real_t * const ANYODE_RESTRICT fy,
                                      bool jok, bool& jac_recomputed, Real_t gamma
                                      )
{
    auto status = AnyODE::Status::success;
    ignore(gamma);
    // See 4.6.9 on page 68 (78) in cvs_guide.pdf (Sundials 2.5)
    if (!jac_cache){
        const int nsat = (geom == Geom::PERIODIC) ? nsidep : 0;
        const int ld = n;
        jac_cache = AnyODE::make_unique<block_diag_ilu::BlockDiagMatrix<Real_t>>(nullptr, N, n, nsidep, nsat, ld);
    }
    if (!jok){
        const int dummy = 0;
        jac_cache->set_to(0);
        status = compressed_jac_cmaj(t, y, fy, jac_cache->m_data, dummy);
        update_prec_cache = true;
        jac_recomputed = true;
    } else jac_recomputed = false;
    nprec_setup++;
    return status;
}
#undef LINC
#undef Y
#undef LAP_WEIGHT
#undef DIV_WEIGHT
#undef GRAD_WEIGHT

template<typename Real_t>
AnyODE::Status
ReactionDiffusion<Real_t>::prec_solve_left(const Real_t t,
                                           const Real_t * const ANYODE_RESTRICT y,
                                           const Real_t * const ANYODE_RESTRICT fy,
                                           const Real_t * const ANYODE_RESTRICT r,
                                           Real_t * const ANYODE_RESTRICT z,
                                           Real_t gamma,
                                           Real_t delta,
                                           const Real_t * const ANYODE_RESTRICT ewt
                                           )
{
    // See 4.6.9 on page 75 in cvs_guide.pdf (Sundials 2.6.2)
    // Solves P*z = r, where P ~= I - gamma*J
    ignore(delta);
    if (ewt)
        throw std::runtime_error("Not implemented.");
    nprec_solve++;

    ignore(t); ignore(fy); ignore(y);
    bool recompute = false;
    if (!prec_cache){
        const int nsat = (geom == Geom::PERIODIC) ? nsidep : 0;
        const int ld = n;
        prec_cache = AnyODE::make_unique<block_diag_ilu::BlockDiagMatrix<Real_t>>(nullptr, N, n, nsidep, nsat, ld);
        recompute = true;
    } else {
        if (update_prec_cache or (old_gamma != gamma))
            recompute = true;
    }
    if (recompute){
        old_gamma = gamma;
        prec_cache->set_to_eye_plus_scaled_mtx(-gamma, *jac_cache);
#if defined(CHEMREAC_WITH_DATA_DUMPING)
        {
            std::ostringstream fname;
            fname << "prec_M_" << std::setfill('0') << std::setw(5) << nprec_solve << ".dat";
            save_array(prec_cache->m_data, prec_cache->m_ndata, fname.str());
        }
#endif
    }

#if defined(CHEMREAC_WITH_DATA_DUMPING)
    {
        std::ostringstream fname;
        fname << "prec_r_" << std::setfill('0') << std::setw(5) << nprec_solve << ".dat";
        save_array(r, n*N, fname.str());
    }
    {
        std::ostringstream fname;
        fname << "prec_g_" << std::setfill('0') << std::setw(5) << nprec_solve << ".dat";
        save_array(&gamma, 1, fname.str());
    }

#endif

    int info;
    if (prec_cache->average_diag_weight(0) > ilu_limit) {
        block_diag_ilu::ILU<Real_t> ilu {*prec_cache};
        nprec_solve_ilu++;
        info = ilu.solve(r, z);
    } else {
        AnyODE::BandedMatrix<Real_t> bm {*prec_cache, get_mlower(), get_mupper()};
        AnyODE::BandedLU<Real_t> lu {&bm};
        lu.factorize();
        nprec_solve_lu++;
        info = lu.solve(r, z);
    }
    if (info == 0)
        return AnyODE::Status::success;
    return AnyODE::Status::recoverable_error;
}

template<typename Real_t>
void
ReactionDiffusion<Real_t>::per_rxn_contrib_to_fi(Real_t t, const Real_t * const ANYODE_RESTRICT y,
                                              int si, Real_t * const ANYODE_RESTRICT out) const
{
    ignore(t);
    auto local_r = AnyODE::buffer_factory<double>(nr);
    fill_local_r_(0, y, AnyODE::buffer_get_raw_ptr(local_r));
    for (int ri=0; ri<nr; ++ri){
	out[ri] = coeff_total[ri*n+si]*local_r[ri];
    }
}

template<typename Real_t>
int
ReactionDiffusion<Real_t>::get_geom_as_int() const
{
    switch(geom){
    case Geom::FLAT :        return 0;
    case Geom::CYLINDRICAL : return 1;
    case Geom::SPHERICAL :   return 2;
    case Geom::PERIODIC :    return 3;
    default:                 return -1;
    }
}

template<typename Real_t>
void
ReactionDiffusion<Real_t>::calc_efield(const Real_t * const linC)
{
    // Prototype for self-generated electric field
    const Real_t F = this->faraday_const; // Faraday's constant
    const Real_t pi = 3.14159265358979324;
    const Real_t eps = eps_rel*vacuum_permittivity;
    Real_t nx, cx = logx ? expb(x[0]) : x[0];
    for (int bi=0; bi<N; ++bi){
        netchg[bi] = 0.0;
        for (int si=0; si<n; ++si)
            netchg[bi] += z_chg[si]*linC[bi*n+si];
    }
    Real_t Q = surf_chg.first;
    for (int bi=0; bi<N; ++bi){
        const Real_t r = logx ? expb(xc[nsidep+bi]) : xc[nsidep+bi];
        nx = logx ? expb(x[bi+1]) : x[bi+1];
        switch(geom){
        case Geom::FLAT:
            efield[bi] = F*Q/eps;
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
        case Geom::PERIODIC:
            efield[bi] = 0;
        }
        cx = nx;
    }
    if (geom == Geom::FLAT){
        Q = surf_chg.second;
        for (int bi=N; bi>0; --bi){
            nx = logx ? expb(x[bi-1]) : x[bi-1];
            efield[bi-1] -= F*Q/eps;
            Q += netchg[bi-1]*(cx - nx);
            cx = nx;
        }
    }
}
} // namespace chemreac

template class chemreac::ReactionDiffusion<double>; // instantiate template
