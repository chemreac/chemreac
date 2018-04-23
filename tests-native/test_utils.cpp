#include "test_utils.h"

using chemreac::ReactionDiffusion;
using std::vector;

// A      -> B
// B + 2C -> A + D

ReactionDiffusion<double> get_four_species_system(int N){
    int n = 4;
    int nr = 2;
    vector<vector<int> > stoich_reac {{0}, {1, 2, 2}};
    vector<vector<int> > stoich_actv;
    vector<vector<int> > stoich_prod {{1}, {1, 3}};
    vector<double> k {0.05, 3.0};
    vector<double> D {.1, .2, .3, .4};
    vector<int> z_chg {0, 0, 0, 0};
    vector<double> mobility {0, 0, 0, 0};
    vector<double> x;
    vector<vector<double> > g_values;
    vector<int> g_value_parents;
    vector<vector<double> > fields;
    vector<int> v;
    int geom = 0;
    bool logy = false, logt = false, logx = false;
    int nstencil = (N == 1) ? 1 : 3;
    for (int ri=0; ri<nr; ++ri)
	stoich_actv.push_back(v);
    for (int i=0; i<=N; ++i)
	x.push_back(1.0 + (double)i*1.0/N);
    return ReactionDiffusion<double>(n, stoich_reac, stoich_prod, k, N, D, z_chg,
                                     mobility, x, stoich_actv, geom, logy, logt, logx, nstencil, true, true);
}
