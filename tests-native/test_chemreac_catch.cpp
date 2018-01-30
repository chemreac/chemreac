#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "catch.hpp"
#include "chemreac/chemreac.hpp"
#include <array>
#include <vector>
#include <utility>

#include "test_utils.h"

#include <iostream>

using std::pair;
using std::vector;
using chemreac::ReactionDiffusion;

TEST_CASE( "jac_times_vec", "[ReactionDiffusion]" ) {

    // this is _get_test_m2 in test_fakelu.py
    auto rd = get_four_species_system(3);

    std::array<double, 3*4> y;
    for (int i=0; i<3; ++i){
        y[4*i + 0] = 1.3;
        y[4*i + 1] = 1e-4;
        y[4*i + 2] = 0.7;
        y[4*i + 3] = 1e-4;
    }

    std::array<double, 3*4> x {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
    std::array<double, 3*4> b;
    std::array<double, 3*4> bref;

    // J*x = b

    std::array<double, 3*4*3*4> J_data;
    auto J = [&](int ri, int ci) -> double& {
        return J_data[ri*3*4+ci];
    };
    rd.dense_jac_rmaj(0.0, &y[0], nullptr, &J_data[0], 3*4);
    for (int ri=0; ri<3*4; ++ri){
        bref[ri] = 0.0;
        for (int ci=0; ci<3*4; ++ci){
            bref[ri] += x[ci]*J(ri, ci);
        }
    }
    rd.jac_times_vec(&x[0], &b[0], 0.0, &y[0], nullptr);
    for (int i=0; i<3*4; ++i)
        REQUIRE( std::abs(bref[i] - b[i]) < 1e-14 );
}

TEST_CASE("private_members", "[ReactionDiffusion]"){
    // from examples/equilibrium.py
    // A + B -> C
    //     C -> A + B
    const int N=1;
    vector<double> D {0.0, 0.0, 0.0};
    vector<int> z_chg {0, 0, 0};
    vector<double> mobility {0.0, 0.0, 0.0};
    vector<double> x {0, 1};
    int geom_=0;
    bool logy=false;
    bool logt=false;
    bool logx=false;
    int nstencil=1;
    double k1=0.9, k2=0.23;
    auto rd = ReactionDiffusion<double>(3, {{0, 1}, {2}}, {{2}, {0, 1}}, {k1, k2}, N, D, z_chg,
                                        mobility, x, {{}, {}}, geom_, logy, logt, logx, nstencil);
    vector<int> ref_ridxs {1, 0};
    vector<int> ref_reaction_orders_seen {1, 2};  // e.g. {1, 3}
    vector<int> ref_n_reac_in_seen_order {1, 1};  // e.g. {2, 1} (that's 3 reactions in total)
    vector<int> ref_active_reac_indices {2, 0, 1}; // e.g. {17, 13, 15, 19, 4} (orders: 1, 1, 3)
    vector<int> ref_n_net_affects {3, 3}; // e.g. {1, 3, 5}
    vector<pair<int, int> > ref_idx_net_stoich {{0, 1}, {1, 1}, {2, -1}, {0, -1}, {1, -1}, {2, 1}};

    REQUIRE( rd.ridxs == ref_ridxs );
    REQUIRE( rd.reaction_orders_seen == ref_reaction_orders_seen );
    REQUIRE( rd.n_reac_in_seen_order == ref_n_reac_in_seen_order );
    REQUIRE( rd.active_reac_indices == ref_active_reac_indices );
    REQUIRE( rd.n_net_affects == ref_n_net_affects );
    REQUIRE( rd.idx_net_stoich == ref_idx_net_stoich );

    double a=3.0, b=5.0, c=7.0;
    vector<double> y {a, b, c};
    vector<double> fout(3);
    rd.rhs(42.0, &y[0], &fout[0]);
    double r1 = k1*y[0]*y[1];
    double r2 = k2*y[2];
    REQUIRE( std::abs(fout[0] - (r2 - r1)) < 1e-14 );
    REQUIRE( std::abs(fout[1] - (r2 - r1)) < 1e-14 );
    REQUIRE( std::abs(fout[2] - (r1 - r2)) < 1e-14 );

    vector<double> jout(9), dfdt(3);
    rd.dense_jac_rmaj(42.0, &y[0], &fout[0], &jout[0], 3, &dfdt[0]);
    vector<double> jref {-k1*b, -k1*a, k2,
            -k1*b, -k1*a, k2,
            k1*b, k1*a, -k2};
    for (int i=0; i<9; ++i)
        REQUIRE( std::abs(jout[i] - jref[i]) < 1e-14 );
}
