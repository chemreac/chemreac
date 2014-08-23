#include <iostream>
#include <vector>
#include <numeric> // std::accumulate
#include <algorithm> // min, max
#include <cassert>
#include "test_utils.h"
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

// A      -> B
// B + 2C -> A + D

int test_f(){
    ReactionDiffusion rd = get_four_species_system(3);
    vector<double> y {1.3, 1e-4, 0.7, 1e-4, 1.3, 1e-4, 0.7, 1e-4, 1.3, 1e-4, 0.7, 1e-4};
    vector<double> ref_f {-0.05*y[0], 0.05*y[0], -2*3.0*y[2]*y[2]*y[1], 3.0*y[2]*y[2]*y[1],\
	    -0.05*y[4], 0.05*y[4], -2*3.0*y[6]*y[6]*y[5], 3.0*y[6]*y[6]*y[5],\
            -0.05*y[8], 0.05*y[8], -2*3.0*y[10]*y[10]*y[9], 3.0*y[10]*y[10]*y[9]};
    double f[12];
    rd.f(0.0, &y[0], f);
    int exit1 = 0;
    for (int i=0; i<12; ++i)
	if (dabs(f[i]-ref_f[i]) > 1e-14){
            std::cout << i << " " << f[i] << " " << ref_f[i] << std::endl;            
            exit1 = 1;
        }
    if (exit1)
        return 1;
    else
        return 0;
}

// n = 4 species
// N = 3 compartments
#define RJ(i, j) ref_jac[(i)*12+j]
int test_dense_jac(){
    ReactionDiffusion rd = get_four_species_system(3);
    const double dx = 1.0 / 3;
    vector<double> y {1.3, 1e-4, 0.7, 1e-4, 1.3, 1e-4, 0.7, 1e-4, 1.3, 1e-4, 0.7, 1e-4};
    double ref_jac[12*12];
    double dense_jac[12*12];

    for (int i=0; i<12*12; ++i){
        ref_jac[i] = 0.0;
        dense_jac[i] = 0.0;
    }
    // First block
    RJ(0,0) = -0.05 -rd.D[0]/(dx*dx);
    RJ(0,4) =  rd.D[0]/(dx*dx);
    RJ(1,0) =  0.05;
    RJ(1,1) =  -rd.D[1]/(dx*dx);
    RJ(1,5) =  rd.D[1]/(dx*dx);
    RJ(2,1) = -2*3.0*y[2]*y[2];
    RJ(2,2) = -2*2*3.0*y[2]*y[1] - rd.D[2]/(dx*dx);
    RJ(2,6) =  rd.D[2]/(dx*dx);
    RJ(3,1) =  3.0*y[2]*y[2];
    RJ(3,2) =  2*3.0*y[2]*y[1];
    RJ(3,3) =  -rd.D[3]/(dx*dx);
    RJ(3,7) =  rd.D[3]/(dx*dx);

    RJ(4,0) = rd.D[0]/(dx*dx);
    RJ(4,4) = -0.05 - 2*rd.D[0]/(dx*dx);
    RJ(4,8) = rd.D[0]/(dx*dx);
    RJ(5,1) = rd.D[1]/(dx*dx);
    RJ(5,4) =  0.05;
    RJ(5,5) = -2*rd.D[1]/(dx*dx);
    RJ(5,9) = rd.D[1]/(dx*dx);
    RJ(6,2) = rd.D[2]/(dx*dx);
    RJ(6,5) = -2*3.0*y[6]*y[6];
    RJ(6,6) = -2*2*3.0*y[6]*y[5]-2*rd.D[2]/(dx*dx);
    RJ(6,10) = rd.D[2]/(dx*dx);
    RJ(7,3) = rd.D[3]/(dx*dx);
    RJ(7,5) =  3.0*y[6]*y[6];
    RJ(7,6) =  2*3.0*y[6]*y[5];
    RJ(7,7) = -2*rd.D[3]/(dx*dx);
    RJ(7,11) = rd.D[3]/(dx*dx);

    RJ(8,4) = rd.D[0]/(dx*dx);
    RJ(8,8) = -0.05 - rd.D[0]/(dx*dx);
    RJ(9,5) = rd.D[1]/(dx*dx);
    RJ(9,8) =  0.05;
    RJ(9,9) = -rd.D[1]/(dx*dx);
    RJ(10,6) = rd.D[2]/(dx*dx);
    RJ(10,9) = -2*3.0*y[10]*y[10];
    RJ(10,10) = -2*2*3.0*y[10]*y[9]-rd.D[2]/(dx*dx);
    RJ(11,7) = rd.D[3]/(dx*dx);
    RJ(11,9) =  3.0*y[10]*y[10];
    RJ(11,10) =  2*3.0*y[10]*y[9];
    RJ(11,11) = -rd.D[3]/(dx*dx);

    rd.dense_jac_rmaj(0.0, &y[0], dense_jac, 12);
    int exit2 = 0;
    for (int i=0; i<12*12; ++i)
	if (dabs(dense_jac[i]-ref_jac[i]) > 1e-14){
            printf("i=%d, dense_jac[i]=%.3f, ref_jac[i]=%.3f\n", i, dense_jac[i], ref_jac[i]);
            exit2 = 1;
        }

    double * bnd_jac = new double[(2*rd.n+1)*(rd.N*rd.n)];
    rd.banded_packed_jac_cmaj(0.0, &y[0], bnd_jac, (2*rd.n+1));
#define BND(i, j) bnd_jac[i-j+rd.n+j*(2*rd.n+1)]
#define DNS(i, j) ref_jac[(i)*rd.n*rd.N+j]
    for (int ri=0; ri<12; ++ri)
	for (int ci=max(0, ri-(int)rd.n); ci<min(rd.n*rd.N, ri+rd.n); ++ci)
	    if (dabs(BND(ri,ci) - DNS(ri,ci)) > 1e-14){
		std::cout << ri << " " << ci << " " << BND(ri,ci) << " " << DNS(ri,ci) << std::endl;
                exit2 = 1;
            }
#undef BND
#undef DNS
    delete []bnd_jac;
    if (exit2)
        return 2;
    else
        return 0;
}
#undef RJ

void bench_f(){
    double t = 0.0;
    int ntimings = 40;
    int N = 500000; // A ridiculous number of bins for 1D but used for benchmarking
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
