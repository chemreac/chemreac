#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "catch.hpp"
#include "chemreac.hpp"
#include <array>

#include "test_utils.h"

#include <iostream>

TEST_CASE( "jac_times_vec", "[ReactionDiffusion]" ) {

    // this is _get_test_m2 in test_fakelu.py
    auto rdp = get_four_species_system(3);
    auto &rd = *rdp;

    std::array<double, 3*4> y;
    for (int i=0; i<3; ++i){
        y[4*i + 0] = 1.3;
        y[4*i + 1] = 1e-4;
        y[4*i + 2] = 0.7;
        y[4*i + 3] = 1e5;
    }

    std::array<double, 3*4> x {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};

    std::array<double, 3*4> f;
    std::array<double, 3*4> fref {
        -0.05*1.3, 0.05*1.3, -2*3.0*1e-4*0.7*0.7, 3.0*1e-4*0.7*0.7,
        -0.05*1.3, 0.05*1.3, -2*3.0*1e-4*0.7*0.7, 3.0*1e-4*0.7*0.7,
        -0.05*1.3, 0.05*1.3, -2*3.0*1e-4*0.7*0.7, 3.0*1e-4*0.7*0.7
    };
    rd.rhs(0.0, y.data(), f.data());
    for (unsigned i=0; i < f.size(); ++i){
        REQUIRE( std::abs(f[i] - fref[i]) < 1e-10 );
    }

    std::array<double, 3*4> b;
    std::array<double, 3*4> bref;

    // J*x = b

    std::array<double, 3*4*3*4> J_data;
    std::memset(J_data.data(), 0, J_data.size()*sizeof(double));
    auto J = [&](int ri, int ci) -> double& {
        return J_data[ri*3*4+ci];
    };
    rd.dense_jac_rmaj(0.0, &y[0], nullptr, &J_data[0], 3*4);
    // Only if D is == 0.0:
    // std::array<double, 3*4*3*4> J_ref {
    //     // import numpy as np
    //     // J_ref = np.zeros((12, 12))
    //     // k1, k2 = 0.05, 3.0
    //     // A, B, C, D = 1.3, 1e-4, 0.7, 1e-4
    //     // j = np.array([[-k1,         0,           0, 0],
    //     //  [k1,          0,           0, 0],
    //     //  [0,   -2*k2*C*C, -2*2*k2*C*B, 0],
    //     //  [0,      k2*C*C,    2*k2*C*B, 0]]
    //     // )
    //     // J_ref[0:4, 0:4] = j
    //     // J_ref[4:8, 4:8] = j
    //     // J_ref[8:12, 8:12] = j
    //     // str(list(J_ref.flat))
    //     -0.05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    //      0.05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    //      0.0, -2.94, -0.00084, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    //      0.0, 1.47, 0.00042, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    //      0.0, 0.0, 0.0, 0.0, -0.05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    //      0.0, 0.0, 0.0, 0.0, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0,
    //      0.0, 0.0, 0.0, 0.0, 0.0, -2.94, -0.00084, 0.0, 0.0, 0.0, 0.0, 0.0,
    //      0.0, 0.0, 0.0, 0.0, 0.0, 1.47, 0.00042, 0.0, 0.0, 0.0, 0.0, 0.0,
    //      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.05, 0.0, 0.0, 0.0,
    //      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05, 0.0, 0.0, 0.0,
    //      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.94, -0.00084, 0.0,
    //      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.47, 0.00042, 0.0
    // };
    // for (unsigned i=0; i < J_ref.size(); ++i){
    //     std::cout << i;
    //     REQUIRE( std::abs( J_ref[i] - J_data[i] ) < 1e-14);
    // }
    // std::array<double, 3*4> bref2 {
    //     -0.1    ,   0.1    ,  -8.82336,   4.41168,  -0.3    ,   0.3    , -20.58672,  10.29336,  -0.5    ,   0.5    , -32.35008,  16.17504
    // };
    for (int ri=0; ri<3*4; ++ri){
        bref[ri] = 0.0;
        for (int ci=0; ci<3*4; ++ci){
            bref[ri] += x[ci]*J(ri, ci);
        }
        // REQUIRE(std::abs(bref[ri] - bref2[ri] < 1e-14));
    }
    rd.jac_times_vec(&x[0], &b[0], 0.0, &y[0], nullptr);
    for (int i=0; i<3*4; ++i){
        std::cout << "jac_times_vec out[i="<< i<<"]=" << b[i] << std::endl;
        REQUIRE( std::abs(bref[i] - b[i]) < 1e-14 );
    }
}
