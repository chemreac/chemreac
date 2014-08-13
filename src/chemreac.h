#ifndef _PVHQOBGMVZECTIJSMOKFUXJXXM
#define _PVHQOBGMVZECTIJSMOKFUXJXXM

#include <vector>
#include <stdexcept>

enum class Geom {FLAT, CYLINDRICAL, SPHERICAL};

namespace chemreac {

using std::vector;

class ReactionDiffusion
{
//private:
public:
    int * coeff_reac;
    int * coeff_prod;
    int * coeff_totl;
    int * coeff_actv;
    double * D_weight; // diffusion weights 
    double * A_weight; // Advection weights 
    vector<uint> i_bin_k;
    uint n_factor_affected_k;
    Geom geom; // Geometry: 0: 1D flat, 1: 1D Cylind, 2: 1D Spherical.

    void _fill_local_r(int, const double * const __restrict__, double * const __restrict__) const;
    void _apply_fd(uint);
    const double * _alloc_and_populate_linC(const double * const __restrict__) const;
    uint _stencil_bi_lbound(uint bi) const;
    uint _xc_bi_map(uint xci) const;

public:
    const uint n; // number of species
    const uint N; // number of compartments
    const uint nstencil; // number of points used in finite difference stencil
    const uint nsidep; // (nstencil-1)/2
    const uint nr; // number of reactions
    const bool logy; // use logarithmic concenctraction
    const bool logt; // use logarithmic time
    const bool logx; // use logarithmic x (space coordinate)
    const vector<vector<uint> > stoich_reac; // Reactants per reaction
    vector<vector<uint> > stoich_actv; // Active reactants per reaction
    const vector<vector<uint> > stoich_prod; // Products per reaction
    vector<double> k; // Rate coefficients (law of mass action)
    vector<double> D; // Diffusion coefficients
    vector<int> z_chg; // ion charge
    vector<double> mobility; // electrical mobility
    const vector<double> x; // Bin edges (length = N+1)
    vector<vector<double> > bin_k_factor; // rate = FACTOR(ri, bi)*k[ri]*C[si1]*C[...]
    vector<uint> bin_k_factor_span; // 
    const bool lrefl, rrefl;
    const bool auto_efield;
    const double surf_chg;
    const double eps;
    double * const efield; // v_d = mu_el*E
    double * xc; // bin centers (length = N+nstencil-1), first bin center: xc[(nstencil-1)/2]

    ReactionDiffusion(uint, 
		      const vector<vector<uint> >, 
		      const vector<vector<uint> >, 
		      vector<double>, 
		      uint,
		      vector<double>, 
                      vector<int>,
                      vector<double>,
		      const vector<double>,
		      vector<vector<uint> >, 
		      vector<vector<double> >,
		      vector<uint>,
		      int geom_=0,
		      bool logy=false,
		      bool logt=false,
                      bool logx=false,
                      uint nstencil=3,
                      bool lrefl=false,
                      bool rrefl=false,
                      bool auto_efield=false,
                      double surf_chg=0.0,
                      double eps=1.0);
    ~ReactionDiffusion();
    void f(double, const double * const, double * const __restrict__);
    void dense_jac_rmaj(double, const double * const, double * const __restrict__, int);
    void dense_jac_cmaj(double, const double * const, double * const __restrict__, int);
    void banded_padded_jac_cmaj(double, const double * const, double * const __restrict__, int);
    void banded_packed_jac_cmaj(double, const double * const, double * const __restrict__, int);
    void per_rxn_contrib_to_fi(double, const double * const __restrict__, uint, double * const __restrict__) const;
    int get_geom_as_int() const;
    void calc_efield(const double * const);
}; // class ReactionDiffusion

} // namespace chemreac
#endif // _PVHQOBGMVZECTIJSMOKFUXJXXM
