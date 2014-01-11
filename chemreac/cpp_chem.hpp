#ifndef _CPP_CHEM_H_
#define _CPP_CHEM_H_

#include <vector>
#include <stdexcept>

using std::vector;

enum class Geom {FLAT, SPHERICAL, CYLINDRICAL};

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
    int mode; // Jacobian storage: 0: DENSE, 1: BANDED, 2: SPARSE
    Geom geom; // Geometry: 0: 1D flat, 1: 1D Spherical, 2: 1D Cylind.
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
		      int);
    ~ReactionDiffusion();
    void f(double, const double * const restrict, double * const restrict) const;
    void dense_jac_rmaj(double, const double * const restrict, double * const restrict, int) const;
    void dense_jac_cmaj(double, const double * const restrict, double * const restrict, int) const;
    void banded_jac_cmaj(double, const double * const restrict, double * const restrict, int) const;
    void banded_packed_jac_cmaj(double, const double * const restrict, double * const restrict, int) const;

};
#endif
