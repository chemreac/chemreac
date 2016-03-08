#include <iostream>
#include <vector>
#include "chemreac.hpp"
#include "cvodes_cxx.hpp"
#include "test_utils.h"

using std::vector;
using chemreac::ReactionDiffusion;

int test_integration(int N){
    ReactionDiffusion<double> rd = get_four_species_system(N);
    vector<double> y;
    for (int i=0; i<N; ++i){
        y.push_back(1.3);
        y.push_back(1e-4);
        y.push_back(0.7);
        y.push_back(1e-4);
    }
    int ny = y.size();
    vector<double> atol {1e-8};
    double rtol {1e-8};
    vector<double> tout {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    double * yout = (double*)malloc(sizeof(double)*tout.size()*ny);
    vector<int> root_indices;
    for (uint i=0; i<tout.size()*ny; ++i) {yout[i]=0.0;}
    cvodes_cxx::simple_predefined<ReactionDiffusion<double> >
        (&rd, atol, rtol, (int)cvodes_cxx::LMM::BDF, &y[0], tout.size(), &tout[0], yout,
         root_indices);
    for (unsigned int tidx=0; tidx<tout.size(); tidx++){
        std::cout << tout[tidx];
        for (int sidx=0; sidx<ny; sidx++){
            std::cout << " " << yout[tidx*ny + sidx];
        }
        std::cout << std::endl;
    }
    free(yout);
    return 0;
}


int main(){
    int status = 0;
    try {
        std::cout << "integrating system...";
        for (int N=1; N<13; N+=2){
            status += test_integration(N);
        }
    } catch (std::exception& e){
        std::cout << e.what() << std::endl;
    }
    return status;
}
