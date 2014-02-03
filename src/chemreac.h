#ifndef _CHEMREAC_H_
#define _CHEMREAC_H_

#include <vector>
#include <stdexcept>

using std::vector;

enum class Geom {FLAT, SPHERICAL, CYLINDRICAL}; // Geom:: -> GEOM_

/* #define Geom int */
/* #define GEOM_FLAT 0 */
/* #define GEOM_SPHERICAL 1 */
/* #define GEOM_CYLINDRICAL 2 */

namespace chemreac {

class ReactionDiffusion
{
private:
    int * coeff_reac;
    int * coeff_prod;
    int * coeff_totl;
    int * coeff_actv;
    double * dx; // bin separations (from center) (length = N)
    vector<int> i_bin_k;
    int n_factor_affected_k;

    void _fill_local_r(int, const double * const restrict, double * const restrict) const;
    double flux(int i, int si, const double * const restrict y) const;
    double diffusion_contrib(int bi, int si, const double * const restrict fluxes) const;
    double diffusion_contrib_jac_prev(int bi) const;
    double diffusion_contrib_jac_next(int bi) const;

public:
    int n; // number of species
    int N; // number of compartments
    int nr; // number of reactions
    Geom geom; // Geometry: 0: 1D flat, 1: 1D Spherical, 2: 1D Cylind.
    int logy; // use logarithmic concenctraction
    int logt; // use logarithmic time
    vector<vector<int> > stoich_reac; // Reactants per reaction
    vector<vector<int> > stoich_actv; // Active reactants per reaction
    vector<vector<int> > stoich_prod; // Products per reaction
    vector<double> k; // Rate coefficients (law of mass action)
    vector<double> D; // Diffusion coefficients
    vector<double> x; // Bin edges (length = N+1)
    vector<vector<double> > bin_k_factor; // rate = FACTOR(ri, bi)*k[ri]*C[si1]*C[...]
    vector<int> bin_k_factor_span; // 

    ReactionDiffusion(int, 
		      vector<vector<int> >, 
		      vector<vector<int> >, 
		      vector<double>, 
		      int,
		      vector<double>, 
		      vector<double>,
		      vector<vector<int> >, 
		      vector<vector<double> >,
		      vector<int>,
		      int,
		      int,
		      int);
    ~ReactionDiffusion();
    void f(double, const double * const restrict, double * const restrict) const;
    void dense_jac_rmaj(double, const double * const restrict, double * const restrict, int) const;
    void dense_jac_cmaj(double, const double * const restrict, double * const restrict, int) const;
    void banded_padded_jac_cmaj(double, const double * const restrict, double * const restrict, int) const;
    void banded_packed_jac_cmaj(double, const double * const restrict, double * const restrict, int) const;

}; // class ReactionDiffusion

}; // namespace chemreac
#endif
