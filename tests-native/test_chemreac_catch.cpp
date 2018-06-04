#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "catch.hpp"
#include "chemreac.hpp"
#include <array>

#include "test_utils.h"

#include <iostream>

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
    for (int i=0; i<3*4; ++i){
        std::cout << "jac_times_vec out[i="<< i<<"]=" << b[i] << std::endl;
        REQUIRE( std::abs(bref[i] - b[i]) < 1e-14 );
    }
}
