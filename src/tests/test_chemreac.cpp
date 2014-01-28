#include <iostream>
#include <vector>
#include <numeric> // std::accumulate
#include <algorithm> // min, max
#include <cassert>
#include "chemreac.h"

#ifdef _OPENMP
  #include <omp.h>
#else
  #include <ctime> // clock_gettime
#endif

double dabs(double a){
    return (a < 0.0) ? -a : a;
}

using std::max;
using std::min;
using std::accumulate;
using std::vector;
using chemreac::ReactionDiffusion;

ReactionDiffusion get_four_species_system(int N){
    int n = 4;
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
    return ReactionDiffusion(\
	n, stoich_reac, stoich_prod, k, N, D, x, stoich_actv,\
	bin_k_factor, bin_k_factor_span, 0);
}

int test_j(){
    ReactionDiffusion rd = get_four_species_system(2);
    vector<double> y0 {1.3, 1e-4, 0.7, 1e-4, 1.3, 1e-4, 0.7, 1e-4};
// about to add test code here
    return 0;
}

int test_f(){
    ReactionDiffusion rd = get_four_species_system(2);
    vector<double> y {1.3, 1e-4, 0.7, 1e-4, 1.3, 1e-4, 0.7, 1e-4};
    vector<double> ref_f {-0.05*y[0], 0.05*y[0], -2*3.0*y[2]*y[2]*y[1], 3.0*y[2]*y[2]*y[1],\
	    -0.05*y[4], 0.05*y[4], -2*3.0*y[6]*y[6]*y[5], 3.0*y[6]*y[6]*y[5]};
    double f[8];
    rd.f(0.0, &y[0], f);
    for (int i=0; i<8; ++i)
	if (dabs(f[i]-ref_f[i]) > 1e-15)
	    return 1;
    return 0;
}

int bench_f(){
    double t = 0.0;
    int ntimings = 10;
    int N = 1000000; // A ridiculous number of bins for 1D but used for benchmarking
    ReactionDiffusion rd = get_four_species_system(N);
    vector<double> y;
    vector<double> b;
    vector<double> timings;
    double timing;
    double best_timing = 1e6;
    double worst_timing = 0.0;
#ifdef _OPENMP
    double t0;
#else
    timespec start, finish;
#endif

    b.reserve(rd.n*N);
    for (auto i = 0; i < N; ++i){
	y.push_back(1.30/(1+1/(i+1)));
	y.push_back(1e-4/(1+1/(i+1)));
	y.push_back(0.70/(1+1/(i+1)));
	y.push_back(1e-4/(1+1/(i+1)));
    }

    for (auto i = 0; i < ntimings; ++i){
	double * const dydt = &b[0];
	const double * const y_ = &y[0];

// Start stop-watch
#ifdef _OPENMP
	t0 = omp_get_wtime();
#else
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
#endif

	rd.f(t, y_, dydt); // heavy lifting

// Stop stop-watch
#ifdef _OPENMP
	timing = omp_get_wtime()-t0;
#else
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &finish);
	timing = (finish.tv_sec-start.tv_sec) + \
	    1e-9*(finish.tv_nsec-start.tv_nsec);
#endif

	best_timing = min(best_timing, timing);
	worst_timing = max(worst_timing, timing);
	timings.push_back(timing);
    }
    std::cout << "Best timing: " << best_timing << std::endl;
    std::cout << "Worst timing: " << worst_timing << std::endl;
    std::cout << "Average timing: " << std::accumulate(timings.begin(), timings.end(), 0.0)/ntimings << std::endl;
    return 0;
}

int main(){
    int status = 0;
    status += test_f();
    status += bench_f();
    status += test_j();
    return status;
}
