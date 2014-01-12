#include <ctime> // clock_gettime
#include <iostream>
#include <vector>
#include <numeric> // std::accumulate
#include <algorithm> // min, max
#include <cassert>
#include "chemreac.h"


using std::max;
using std::min;
using std::accumulate;
using std::vector;


void test_f(){
    double t = 0.0;
    int ntimings = 10;
    int n = 4;
    int N = 400000; // A ridiculous number of bins for 1D but used for benchmarking
    vector<vector<int> > stoich_reac {{0}, {1, 2, 2}};
    vector<vector<int> > stoich_actv;
    vector<vector<int> > stoich_prod {{1}, {1, 3}};
    vector<double> k {0.05, 3.0};
    vector<double> D {.1, .2, .3, .4};
    vector<double> x;
    vector<vector<double> > bin_k_factor;
    vector<int> bin_k_factor_span;
    vector<int> v;
    for (int i=0; i<2; ++i)
	stoich_actv.push_back(v);
    for (int i=0; i<N+1; ++i)
	x.push_back((double)i);
    printf("x.size()=%d\n", x.size());
    chemreac::ReactionDiffusion rd(n, stoich_reac, stoich_prod, k, N, D, x, stoich_actv,\
				   bin_k_factor, bin_k_factor_span, 0, 0);
    vector<double> y;
    vector<double> b;
    vector<double> timings;
    double timing;
    double best_timing = 1e6;
    double worst_timing = 0.0;
    timespec start, finish;

    b.reserve(n*N);
    for (auto i = 0; i < N; ++i){
	y.push_back(1.30/(1+1/(i+1)));
	y.push_back(1e-4/(1+1/(i+1)));
	y.push_back(0.70/(1+1/(i+1)));
	y.push_back(1e-4/(1+1/(i+1)));
    }

    for (auto i = 0; i < ntimings; ++i){
	double * const dydt = &b[0];
	const double * const y_ = &y[0];
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
	rd.f(t, y_, dydt); // heavy lifting
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &finish);
	timing = (finish.tv_sec-start.tv_sec) + \
	    1e-9*(finish.tv_nsec-start.tv_nsec);
	best_timing = min(best_timing, timing);
	worst_timing = max(worst_timing, timing);
	timings.push_back(timing);
    }
    std::cout << "Best timing: " << best_timing << std::endl;
    std::cout << "Worst timing: " << worst_timing << std::endl;
    std::cout << "Average timing: " << std::accumulate(timings.begin(), timings.end(), 0.0)/ntimings << std::endl;

}

int main(){
    test_f();
    return 0;
}
