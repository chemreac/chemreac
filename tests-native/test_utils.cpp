#include "test_utils.h"

using chemreac::ReactionDiffusion;
using std::vector;

// A      -> B       k0
// B + 2C -> B + D   k1

// f:
// [-k0*A, k0*A, -2*k1*B*C*C, k1*B*C*C]

// J:
// -k0   0          0          0
// k0    0          0          0
//  0   -2*C*C*k1  -4*C*B*k1   0
//  0    k1*C*C     2*C*B*k1   0


std::unique_ptr<ReactionDiffusion<double>> get_four_species_system(int N){
    int n = 4;
    int nr = 2;
    vector<vector<int> > stoich_actv {{0}, {1, 2, 2}};
    vector<vector<int> > stoich_inact;
    vector<vector<int> > stoich_prod {{1}, {1, 3}};
    vector<double> k {0.05, 3.0};
    vector<double> D(N*n);
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
    for (int bi=0; bi<N; ++bi){
        D[bi*n + 0] = .1;
        D[bi*n + 1] = .2;
        D[bi*n + 2] = .3;
        D[bi*n + 3] = .4;
    }
    for (int ri=0; ri<nr; ++ri)
	stoich_inact.push_back(v);
    for (int i=0; i<=N; ++i)
	x.push_back(1.0 + (double)i*1.0/N);
    return AnyODE::make_unique<ReactionDiffusion<double>>(
        n, stoich_actv, stoich_prod, k, N, D, z_chg,
        mobility, x, stoich_inact, geom, logy, logt, logx, nstencil, true, true);
}
