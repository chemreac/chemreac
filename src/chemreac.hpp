#ifndef CHEMREAC_PVHQOBGMVZECTIJSMOKFUXJXXM
#define CHEMREAC_PVHQOBGMVZECTIJSMOKFUXJXXM

#include <vector>
#include <utility>
#include <stdexcept>
#include <memory> // unique_ptr
#include "block_diag_ilu.hpp"

enum class Geom {FLAT, CYLINDRICAL, SPHERICAL};

namespace chemreac {

using std::vector;
using std::pair;
using block_diag_ilu::make_unique;

class ReactionDiffusion
{
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
    const pair<double, double> surf_chg;
    const double eps;
    const double faraday_const;
    double * const efield; // v_d = mu_el*E
    double * const netchg;
private:
    block_diag_ilu::BlockDiagMat *jac_cache {nullptr};
    block_diag_ilu::BlockDiagMat *prec_cache {nullptr};

public:
    double * xc; // bin centers (length = N+nstencil-1), first bin center: xc[(nstencil-1)/2]
    long neval_f {0};
    long neval_j {0};
    long nprec_setup {0};
    long nprec_solve {0};
    long njacvec_dot {0};

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
                      pair<double, double> surf_chg={0, 0},
                      double eps=1.0,
                      double faraday_const=96485.3399); //  [C/mol]
    ~ReactionDiffusion();
    void f(double, const double * const, double * const __restrict__);
    void dense_jac_rmaj(double, const double * const __restrict__, const double * const __restrict__, double * const __restrict__, int);
    void dense_jac_cmaj(double, const double * const __restrict__, const double * const __restrict__, double * const __restrict__, int);
    void banded_padded_jac_cmaj(double, const double * const __restrict__, const double * const __restrict__, double * const __restrict__, int);
    void banded_packed_jac_cmaj(double, const double * const __restrict__,  const double * const __restrict__, double * const __restrict__, int);
    void compressed_jac_cmaj(double, const double * const __restrict__, const double * const __restrict__, double * const __restrict__, int);
    // For iterative linear solver
    // void local_reaction_jac(const uint, const double * const, double * const __restrict__, double) const;
    void jac_times_vec(const double * const __restrict__ vec,
                       double * const __restrict__ out,
                       double t, const double * const __restrict__ y,
                       const double * const __restrict__ fy
                       );
    void prec_setup(double t, const double * const __restrict__ y, 
                    const double * const __restrict__ fy, 
                    bool jok, bool& jac_recomputed, double gamma);
    void prec_solve_left(const double t, const double * const __restrict__ y,
                         const double * const __restrict__ fy, 
                         const double * const __restrict__ r, 
                         double * const __restrict__ z,
                         double gamma);

    void per_rxn_contrib_to_fi(double, const double * const __restrict__, uint, double * const __restrict__) const;
    int get_geom_as_int() const;
    void calc_efield(const double * const);
}; // class ReactionDiffusion

} // namespace chemreac
#endif // CHEMREAC_PVHQOBGMVZECTIJSMOKFUXJXXM
