#include <iostream>
#include <vector>
#include <memory>
#include <numeric> // std::accumulate
#include <algorithm> // min, max
#include <cassert>
#include "chemreac.hpp"
#include "test_utils.h"

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

int test_rhs(){
    auto rd = get_four_species_system(3);
    vector<double> y {1.3, 1e-4, 0.7, 1e-4, 1.3, 1e-4, 0.7, 1e-4, 1.3, 1e-4, 0.7, 1e-4};
    vector<double> ref_f {-0.05*y[0], 0.05*y[0], -2*3.0*y[2]*y[2]*y[1], 3.0*y[2]*y[2]*y[1],\
	    -0.05*y[4], 0.05*y[4], -2*3.0*y[6]*y[6]*y[5], 3.0*y[6]*y[6]*y[5],\
            -0.05*y[8], 0.05*y[8], -2*3.0*y[10]*y[10]*y[9], 3.0*y[10]*y[10]*y[9]};
    double f[12];
    rd->rhs(0.0, &y[0], f);
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
int test_jac(){
    auto rdp = get_four_species_system(3);
    auto &rd = *rdp;

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

    rd.dense_jac_rmaj(0.0, &y[0], nullptr, dense_jac, 12);
    int exit2 = 0;
    for (int i=0; i<12*12; ++i)
	if (dabs(dense_jac[i]-ref_jac[i]) > 1e-14){
            printf("i=%d, dense_jac[i]=%.3f, ref_jac[i]=%.3f\n", i, dense_jac[i], ref_jac[i]);
            exit2 = exit2 | 1;
        }

    // Banded jacobian
    const int ld = (3*rd.n_jac_diags*rd.n + 1);
    const int ny = rd.N*rd.n;
    double * bnd_jac = new double[ld*(rd.N*rd.n)];
    for (int i=0; i<ld*ny; ++i) bnd_jac[i] = 0.0;
    rd.banded_jac_cmaj(0.0, &y[0], nullptr, bnd_jac + rd.n*rd.n_jac_diags, ld);
#define BND(i, j) bnd_jac[i-j+2*rd.n+j*ld]
#define DNS(i, j) ref_jac[(i)*rd.n*rd.N+j]
    for (int ri=0; ri<12; ++ri)
	for (int ci=std::max(0, ri-(int)rd.n); ci<std::min(rd.n*rd.N, ri+rd.n); ++ci)
	    if (dabs(BND(ri,ci) - DNS(ri,ci)) > 1e-14){
		std::cout << ri << " " << ci << " " << BND(ri,ci) << " " << DNS(ri,ci) << std::endl;
                exit2 = exit2 | 2;
            }
#undef BND
    delete []bnd_jac;

    // Compressed jacobian
    vector<double> cmprs_jac(rd.n*rd.n*rd.N + 2*rd.n*(rd.N-1), 0);
    rd.compressed_jac_cmaj(0.0, &y[0], nullptr, &cmprs_jac[0], rd.n);
    // std::cout << "n_jac_diags = " << rd.n_jac_diags << std::endl;
#define CMPRS(bi, ri, ci) cmprs_jac[bi*rd.n*rd.n + ci*rd.n + ri]
#define SUB(bi, ci) cmprs_jac[rd.N*rd.n*rd.n + rd.n*bi + ci]
#define SUP(bi, ci) cmprs_jac[rd.N*rd.n*rd.n + (rd.N-1)*rd.n + rd.n*bi + ci]
    // diagonal blocks
    for (int bi=0; bi<rd.N; ++bi)
        for (int ci=0; ci<rd.n; ++ci)
            for (int ri=0; ri<rd.n; ++ri)
                if (dabs(CMPRS(bi, ri, ci) - DNS(bi*rd.n + ri, bi*rd.n + ci)) > 1e-14){
                    std::cout << "CMPRS: " << bi << " " << ci << " " << ri << " " <<
                        CMPRS(bi, ri, ci) << " " << DNS(bi*rd.n + ri, bi*rd.n + ci) << std::endl;
                    exit2 = exit2 | 4;
                }
    for (int bi=0; bi<rd.N-1; ++bi)
        for (int ci=0; ci<rd.n; ++ci){
            // sub diagonal
            if (dabs(SUB(bi, ci) - DNS((bi+1)*rd.n + ci, bi*rd.n + ci)) > 1e-14){
                std::cout << "SUB: " << bi << " " << ci << " " << SUB(bi, ci) << " " << DNS((bi+1)*rd.n + ci, bi*rd.n + ci) << std::endl;
                exit2 = exit2 | 4;
            }
            // sup diagonal
            if (dabs(SUP(bi, ci) - DNS(bi*rd.n + ci, (bi+1)*rd.n + ci)) > 1e-14){
                std::cout << "SUP:" << bi << " " << ci << " " << SUB(bi, ci) << " " << DNS((bi+1)*rd.n + ci, bi*rd.n + ci) << std::endl;
                exit2 = exit2 | 4;
            }
        }
#undef SUP
#undef SUP
#undef CMPRS

    vector<double> vec(12);
    for (int i=0; i<12; ++i) vec[i] = i+2.0;
    vector<double> out(12);
    for (int i=0; i<12; ++i) out[i] = 0.0;
    rd.jac_times_vec(&vec[0], &out[0], 0.0, &y[0], nullptr);
    for (int ri=0; ri<12; ++ri){
        double val = 0.0;
        for (int ci=0; ci<12; ++ci){
            // std::cout << "ri=" << ri << " ci=" << ci << " => add: " <<
            //     DNS(ri, ci)*vec[ci] << std::endl;
            val += DNS(ri, ci)*vec[ci];
        }
        if (dabs(out[ri]-val) > 1e-13){
            std::cout << "jac_times_vec failed for ri=" << ri << " ref[ri]=" << val <<
                " out[ri]=" << out[ri] << " (difference = " << out[ri]-val << ")" << std::endl;
            exit2 = exit2 | 8;
        } else {
            std::cout << "jac_times_vec   ok   for ri=" << ri << " ref[ri]=" << val <<
                " out[ri]=" << out[ri] << " (difference = " << out[ri]-val << ")" << std::endl;
        }

    }
#undef DNS

    if (exit2){
        printf("exit2: %d\n", exit2);
        return 1;
    }
    return 0;

}
#undef RJ

void bench_rhs(){
    double t = 0.0;
    int ntimings = 40;
    int N = 500000; // A ridiculous number of bins for 1D but used for benchmarking
    auto rdp = get_four_species_system(N);
    auto &rd = *rdp;
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

	rd.rhs(t, y_, dydt); // heavy lifting

// Stop stop-watch
#ifdef _OPENMP
	timing = omp_get_wtime()-t0;
#else
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &finish);
	timing = (finish.tv_sec-start.tv_sec) + \
	    1e-9*(finish.tv_nsec-start.tv_nsec);
#endif

	best_timing = std::min(best_timing, timing);
	worst_timing = std::max(worst_timing, timing);
	timings.push_back(timing);
    }
    std::cout << "Best timing: " << best_timing << std::endl;
    std::cout << "Worst timing: " << worst_timing << std::endl;
    std::cout << "Average timing: " << std::accumulate(timings.begin(), timings.end(), 0.0)/ntimings << std::endl;
}

std::unique_ptr<ReactionDiffusion<double>> _get_single_specie_system(int N, int z){
    int n = 1;
    vector<vector<int> > stoich_reac {};
    vector<vector<int> > stoich_actv {};
    vector<vector<int> > stoich_prod {};
    vector<double> k {};
    vector<double> D(N, 1.0);
    vector<int> z_chg {z};
    vector<double> mobility {1.0};
    vector<double> x;
    vector<int> v;
    int geom = 0;
    bool logy = false, logt = false, logx = false;
    int nstencil = (N == 1) ? 1 : 3;
    for (int i=0; i<=N; ++i)
	x.push_back(1.0 + (double)i*1.0/N);
    std::pair<double, double> surf_chg(0, 0);
    return AnyODE::make_unique<ReactionDiffusion<double>>(
        n, stoich_reac, stoich_prod, k, N, D, z_chg,
        mobility, x, stoich_actv, geom, logy, logt, logx,
        nstencil, true, true, true, surf_chg, 1.0);
}

int test_calc_efield(){
    auto rdp = _get_single_specie_system(5, 1);
    auto &rd = *rdp;
    vector<double> y {1.0, 2.0, 3.0, 2.0, 1.0};
    const double factor = 0.2*96485.3399/8.854187817e-12;
    vector<double> ref_efield {-8*factor, -5*factor, 0, 5*factor, 8*factor};
    rd.calc_efield(&y[0]);
    int fail = 0;
    for (unsigned i=0; i<ref_efield.size(); ++i)
	if (dabs((rd.efield[i]-ref_efield[i])/ref_efield[i]) > 1e-10){
            std::cout << i << " " << rd.efield[i] << " " << ref_efield[i] << std::endl;
            fail = 1;
        }
    if (fail)
        return 4;
    else
        return 0;
}

int main(){
    int status = 0;
    try {
        std::cout << "test_f..." << std::endl;
        status += test_rhs();
        std::cout << "test_jac..." << std::endl;
        status += test_jac();
        std::cout << "test_calc_efield..."  << std::endl;
        status += test_calc_efield();
#ifdef BENCHMARK
        std::cout << "bench_f..." << std::endl;
        bench_rhs();
#endif
    } catch (std::exception& e){
        std::cout << e.what() << std::endl;
    }
    return status;
}
