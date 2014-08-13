#include <iostream>
#include <vector>
#include "chemreac.h"
#include "chemreac_sundials.h"
#include "test_utils.h"

using std::vector;

int test_integration(){
    ReactionDiffusion rd = get_four_species_system(1);
    vector<double> y {1.3, 1e-4, 0.7, 1e-4};
    int ny = y.size();
    vector<double> atol {1e-8, 1e-8, 1e-8, 1e-8};
    double rtol {1e-8};
    vector<double> tout {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    double * yout = (double*)malloc(sizeof(double)*tout.size()*ny);
    for (int i=0; i<tout.size()*ny; ++i) {yout[i]=0.0;}
    chemreac_sundials::direct<double, ReactionDiffusion>(&rd, atol, rtol, 1, 
                                                         &y[0], tout.size(), &tout[0], yout);
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
        status += test_integration();
    } catch (std::exception& e){
        std::cout << e.what() << std::endl;
    }
    return status;
}
