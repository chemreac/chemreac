#ifndef _CHEMREAC_H_
#define _CHEMREAC_H_

#include <vector>
#include <stdexcept>

using std::vector;

enum class Geom {FLAT, CYLINDRICAL, SPHERICAL}; // Geom:: -> GEOM_

namespace chemreac {

#define ALIGNED16 __attribute__((aligned(16)))


class ReactionDiffusion
{
//private:
public:
    int * coeff_reac;
    int * coeff_prod;
    int * coeff_totl;
    int * coeff_actv;
    double * D_weight; // diffusion weights 
    vector<uint> i_bin_k;
    uint n_factor_affected_k;
    Geom geom; // Geometry: 0: 1D flat, 1: 1D Cylind, 2: 1D Spherical.

    void _fill_local_r(int, const double * const restrict, double * const restrict) const;
    double flux(int i, int si, const double * const restrict y) const;
    double diffusion_contrib(int bi, int si, const double * const restrict fluxes) const;
    double diffusion_contrib_jac_prev(int bi) const;
    double diffusion_contrib_jac_next(int bi) const;
    void _apply_fd(int, int);
    double * const restrict liny ALIGNED16;
    const int nliny;
    
public:
    const uint n; // number of species
    const uint N; // number of compartments
    const uint nstencil; // number of points used in finite difference stencil
    const uint nr; // number of reactions
    bool logy; // use logarithmic concenctraction
    bool logt; // use logarithmic time
    const vector<vector<uint> > stoich_reac; // Reactants per reaction
    vector<vector<uint> > stoich_actv; // Active reactants per reaction
    const vector<vector<uint> > stoich_prod; // Products per reaction
    vector<double> k; // Rate coefficients (law of mass action)
    vector<double> D; // Diffusion coefficients
    const vector<double> x; // Bin edges (length = N+1)
    vector<vector<double> > bin_k_factor; // rate = FACTOR(ri, bi)*k[ri]*C[si1]*C[...]
    vector<uint> bin_k_factor_span; // 
    double * xc; // bin centers (length = N-1)

    ReactionDiffusion(uint, 
		      const vector<vector<uint> >, 
		      const vector<vector<uint> >, 
		      vector<double>, 
		      uint,
		      vector<double>, 
		      const vector<double>,
		      vector<vector<uint> >, 
		      vector<vector<double> >,
		      vector<uint>,
		      int,
		      bool,
		      bool,
                      uint nstencil = 3);
    ~ReactionDiffusion();
    void f(double, const double * const restrict, double * const restrict) const;
    void dense_jac_rmaj(double, const double * const restrict, double * const restrict, int) const;
    void dense_jac_cmaj(double, const double * const restrict, double * const restrict, int) const;
    void banded_padded_jac_cmaj(double, const double * const restrict, double * const restrict, int) const;
    void banded_packed_jac_cmaj(double, const double * const restrict, double * const restrict, int) const;
    void per_rxn_contrib_to_fi(double, const double * const restrict, uint, double * const restrict) const;
    int get_geom_as_int() const;

}; // class ReactionDiffusion

}; // namespace chemreac
#endif
