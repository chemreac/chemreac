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
    int nr = 2;
    vector<vector<int> > stoich_reac {{0}, {1, 2, 2}};
    vector<vector<int> > stoich_actv;
    vector<vector<int> > stoich_prod {{1}, {1, 3}};
    vector<double> k {0.05, 3.0};
    vector<double> D {.1, .2, .3, .4};
    vector<double> x;
    vector<vector<double> > bin_k_factor;
    vector<int> bin_k_factor_span;
    vector<int> v;
    for (int ri=0; ri<nr; ++ri)
	stoich_actv.push_back(v);
    for (int i=0; i<N+1; ++i)
	x.push_back((double)i);
    return ReactionDiffusion(\
	n, stoich_reac, stoich_prod, k, N, D, x, stoich_actv,\
	bin_k_factor, bin_k_factor_span, 0, 0, 0);
}

#define RJ(i, j) ref_jac[(i)*8+j]
int test_dense_jac(){
    ReactionDiffusion rd = get_four_species_system(2);
    vector<double> y {1.3, 1e-4, 0.7, 1e-4, 1.3, 1e-4, 0.7, 1e-4};
    double ref_jac[8*8];
    for (int i=0; i<8*8; ++i)
        ref_jac[i] = 0.0;
    // First block
    RJ(0,0) = -0.05 -rd.D[0];
    RJ(0,4) =  rd.D[0];
    RJ(1,0) =  0.05;
    RJ(1,1) =  -rd.D[1];
    RJ(1,5) =  rd.D[1];
    RJ(2,1) = -2*3.0*y[2]*y[2];
    RJ(2,2) = -2*2*3.0*y[2]*y[1] - rd.D[2];
    RJ(2,6) =  rd.D[2];
    RJ(3,1) =  3.0*y[2]*y[2];
    RJ(3,2) =  2*3.0*y[2]*y[1];
    RJ(3,3) =  -rd.D[3];
    RJ(3,7) =  rd.D[3];

    RJ(4,0) = rd.D[0];
    RJ(4,4) = -0.05 - rd.D[0];
    RJ(5,1) = rd.D[1];
    RJ(5,4) =  0.05;
    RJ(5,5) = -rd.D[1];
    RJ(6,2) = rd.D[2];
    RJ(6,5) = -2*3.0*y[6]*y[6];
    RJ(6,6) = -2*2*3.0*y[6]*y[5]-rd.D[2];
    RJ(7,3) = rd.D[3];
    RJ(7,5) =  3.0*y[6]*y[6];
    RJ(7,6) =  2*3.0*y[6]*y[5];
    RJ(7,7) = -rd.D[3];

    double dense_jac[8*8];
    rd.dense_jac_rmaj(0.0, &y[0], dense_jac, 8);
    for (int i=0; i<8*8; ++i)
	if (dabs(dense_jac[i]-ref_jac[i]) > 1e-15)
	    return 2;

    double * bnd_jac = new double[(2*rd.n+1)*(rd.N*rd.n)];
    rd.banded_packed_jac_cmaj(0.0, &y[0], bnd_jac, (2*rd.n+1));
#define BND(i, j) bnd_jac[i-j+rd.n+j*(2*rd.n+1)]
#define DNS(i, j) ref_jac[(i)*rd.n*rd.N+j]
    for (int ri=0; ri<8; ++ri)
	for (int ci=max(0, ri-rd.n); ci<min(rd.n*rd.N, ri+rd.n); ++ci)
	    if (dabs(BND(ri,ci) - DNS(ri,ci)) > 1e-15)
		std::cout << ri << " " << ci << " " << BND(ri,ci) << " " << DNS(ri,ci) << std::endl;//return 1;
#undef BND
#undef DNS
    delete []bnd_jac;
    return 0;
}
#undef RJ

int test_f(){
    ReactionDiffusion rd = get_four_species_system(3);
    vector<double> y {1.3, 1e-4, 0.7, 1e-4, 1.3, 1e-4, 0.7, 1e-4, 1.3, 1e-4, 0.7, 1e-4};
    vector<double> ref_f {-0.05*y[0], 0.05*y[0], -2*3.0*y[2]*y[2]*y[1], 3.0*y[2]*y[2]*y[1],\
	    -0.05*y[4], 0.05*y[4], -2*3.0*y[6]*y[6]*y[5], 3.0*y[6]*y[6]*y[5],\
            -0.05*y[8], 0.05*y[8], -2*3.0*y[10]*y[10]*y[9], 3.0*y[10]*y[10]*y[9]};
    double f[12];
    rd.f(0.0, &y[0], f);
    for (int i=0; i<12; ++i)
	if (dabs(f[i]-ref_f[i]) > 1e-15)
	    return 1;
    return 0;
}

void bench_f(){
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
}

int main(){
    int status = 0;
    try {
        std::cout << "test_f...";
        status += test_f();
        std::cout << std::endl <<"test_dense_jac...";
        status += test_dense_jac();
        std::cout << std::endl << "bench_f...";
        bench_f();
    } catch (std::exception& e){
        std::cout << e.what() << std::endl;
    }
    return status;
}
