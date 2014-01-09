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
    vector<double> dx; // bin separations (from center) (length = N)

    void _fill_local_r(double*, double*);
    double flux(int bi, int si, const double * const restrict y);
    double diffusion_contrib(int bi, int si, const double * const restrict fluxes);
    double diffusion_contrib_jac_prev(int bi, int si, const double * const restrict y);
    double diffusion_contrib_jac_next(int bi, int si, const double * const restrict y);

public:
    int n; // number of species
    int N; // number of compartments
    int nr; // number of reactions
    int mode; // Jacobian storage: 0: DENSE, 1: BANDED, 2: SPARSE
    Geom geom; // Geometry: 0: 1D flat, 1: 1D Spherical, 2: 1D Cylind.
    int nfeval; // Number of funcion evaluations
    int njeval; // Number of Jacobian evaluations
    vector<vector<int> > stoich_reac; // Reactants per reaction
    vector<vector<int> > stoich_actv; // Active reactants per reaction
    vector<vector<int> > stoich_prod; // Products per reaction
    vector<double> k; // Rate coefficients (law of mass action)
    vector<double> D; // Diffusion coefficients
    vector<double> x; // Bin edges (length = N+1)

    ReactionDiffusion(int, 
		      int,
		      vector<vector<int> >, 
		      vector<vector<int> >, 
		      vector<vector<int> >, 
		      vector<double>, 
		      vector<double>, 
		      vector<double>,
		      int,
		      int);
    ~ReactionDiffusion();
    void f(double, double*, double*);
    void dense_jac_rmaj(double, double*, double*, int);
    void dense_jac_cmaj(double, double*, double*, int);
    void banded_jac_cmaj(double, double*, double*, int);
    void banded_packed_jac_cmaj(double, double*, double*, int);

};
#endif
