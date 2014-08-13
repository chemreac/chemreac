#include "test_utils.h"

using chemreac::ReactionDiffusion;
using std::vector;

ReactionDiffusion get_four_species_system(int N){
    uint n = 4;
    uint nr = 2;
    vector<vector<uint> > stoich_reac {{0}, {1, 2, 2}};
    vector<vector<uint> > stoich_actv;
    vector<vector<uint> > stoich_prod {{1}, {1, 3}};
    vector<double> k {0.05, 3.0};
    vector<double> D {.1, .2, .3, .4};
    vector<int> z_chg {0, 0, 0, 0};
    vector<double> mobility {0, 0, 0, 0};
    vector<double> x;
    vector<vector<double> > bin_k_factor;
    vector<uint> bin_k_factor_span;
    vector<uint> v;
    for (int ri=0; ri<nr; ++ri)
	stoich_actv.push_back(v);
    for (int i=0; i<N+1; ++i)
	x.push_back((double)i);
    return ReactionDiffusion(n, stoich_reac, stoich_prod, k, N, D, z_chg,
                             mobility, x, stoich_actv,	bin_k_factor,
                             bin_k_factor_span, 0, false, false, false, 1, true, true);
}
