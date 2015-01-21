#ifndef BLOCK_DIAG_ILU_GOB3CSYR2HBHUEX4HJGA3584
#define BLOCK_DIAG_ILU_GOB3CSYR2HBHUEX4HJGA3584
#include <type_traits>
#include <utility>
#include <memory>
#include <cmath> // std::abs

// block_diag_ilu
// ==============
// Algorithm: Incomplete LU factorization of block diagonal matrices with weak sub-/super-diagonals
// Language: C++11
// License: Open Source, see LICENSE.txt (BSD 2-Clause license)
// Author: Björn Dahlgren 2015
// URL: https://github.com/chemreac/block_diag_ilu


namespace block_diag_ilu {

    // make_unique<T[]>() only in C++14, work around:
    // begin copy paste from http://stackoverflow.com/a/10150181/790973
    template <class T, class ...Args>
    typename std::enable_if
    <
        !std::is_array<T>::value,
        std::unique_ptr<T>
        >::type
    make_unique(Args&& ...args)
    {
        return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
    }

    template <class T>
    typename std::enable_if
    <
        std::is_array<T>::value,
        std::unique_ptr<T>
        >::type
    make_unique(std::size_t n)
    {
        typedef typename std::remove_extent<T>::type RT;
        return std::unique_ptr<T>(new RT[n]);
    }
    // end copy paste from http://stackoverflow.com/a/10150181/790973

#if defined(WITH_BLOCK_DIAG_ILU_DGETRF)
    inline int dgetrf_square(const int dim, double * const __restrict__ a,
                      const int lda, int * const __restrict__ ipiv);
#else
    extern "C" void dgetrf_(const int* dim1, const int* dim2, double* a, int* lda, int* ipiv, int* info);
#endif
    inline void rowpiv2rowbycol(int n, const int * const piv, int * const rowbycol) {
        for (int i = 0; i < n; ++i){
            rowbycol[i] = i;
        }
        for (int i=0; i<n; ++i){
            int j = piv[i] - 1; // Fortran indexing starts with 1
            if (i != j){
                int tmp = rowbycol[j];
                rowbycol[j] = rowbycol[i];
                rowbycol[i] = tmp;
            }
        }
    }

    inline void rowbycol2colbyrow(int n, const int * const rowbycol, int * const colbyrow){
        for (int i=0; i<n; ++i){
            for (int j=0; j<n; ++j){
                if (rowbycol[j] == i){
                    colbyrow[i] = j;
                    break;
                }
            }
        }
    }

    constexpr int diag_store_len(int N, int n, int ndiag) {
        return n*(N*ndiag - (ndiag*ndiag + ndiag)/2);
    }

    class ILU {
        double * const __restrict__ block_data;
        double * const __restrict__ sub_data;
        const double * const __restrict__ sup_data;
    public:
        const int nblocks, blockw, ndiag;
    private:
        int ld_block_data;
        std::unique_ptr<int[]> piv, rowbycol, colbyrow;

        inline double lu_get(const int blocki, const int rowi, 
                             const int coli) const {
            const int blockspan = (this->ld_block_data)*(this->blockw);
            return this->block_data[blocki*blockspan + coli*(this->ld_block_data) + rowi];
        }
    public:
        inline double sub_get(const int diagi, const int blocki,
                              const int coli) const {
            const int skip_ahead = diag_store_len(this->nblocks, this->blockw, diagi);
            return this->sub_data[skip_ahead + blocki*(this->blockw) + coli];
        }
        inline double sup_get(const int diagi, const int blocki, const int coli) const {
            const int skip_ahead = diag_store_len(this->nblocks, this->blockw, diagi);
            return this->sup_data[skip_ahead + blocki*(this->blockw) + coli];
        }
#ifdef UNIT_TEST
        int piv_get(const int idx) { return this->piv[idx]; }
        int rowbycol_get(const int idx) { return this->rowbycol[idx]; }
        int colbyrow_get(const int idx) { return this->colbyrow[idx]; }
#endif
        ILU(double * const __restrict__ block_data,
            double * const __restrict__ sub_data,
            const double * const __restrict__ sup_data,
            int nblocks, 
            int blockw, int ndiag,
            int ld_block_data=0)
            // block_data : column major ordering (Fortran style), 
            //     (will be overwritten)
            // sub_data : sub[0], sub[1], ... sup[ndiag-1]
            // sup_data : sup[0], sup[1], ... sup[ndiag-1]
            // use ld_block_data to avoid false sharing in parallelized
            // execution (TODO: guard against false sharing in diagdata..)
            : block_data(block_data), sub_data(sub_data), sup_data(sup_data),
              nblocks(nblocks), blockw(blockw), ndiag(ndiag), 
              ld_block_data((ld_block_data) ? ld_block_data : blockw),
              piv(make_unique<int[]>(blockw*nblocks)),
              rowbycol(make_unique<int[]>(blockw*nblocks)),
              colbyrow(make_unique<int[]>(blockw*nblocks)) {
            int info_ = 0;
#if defined(WITH_BLOCK_DIAG_ILU_OPENMP)
#pragma omp parallel for
#endif
            for (int bi=0; bi<nblocks; ++bi){
#if defined(WITH_BLOCK_DIAG_ILU_DGETRF)
                int info = dgetrf_square(
                           this->blockw,
                           &block_data[bi*blockw*(this->ld_block_data)],
                           this->ld_block_data,
                           &(this->piv[bi*blockw]));
#else
                int info;
                dgetrf_(&(this->blockw),
                        &(this->blockw),
                        &block_data[bi*blockw*(this->ld_block_data)],
                        &(this->ld_block_data),
                        &(this->piv[bi*blockw]),
                        &info);
#endif
                if ((info != 0) && (info_ == 0))
                    info_ = info;
                for (int ci = 0; ci < blockw; ++ci){
                    for (int di = 0; (di < (this->ndiag)) && (bi+di < (this->nblocks) - 1); ++di){
                        const int skip_ahead = diag_store_len(this->nblocks, this->blockw, di);
                        const int gi = skip_ahead + bi*(this->blockw) + ci;
                        this->sub_data[gi] = sub_data[gi]/(this->lu_get(bi, ci, ci));
                    }
                }
                rowpiv2rowbycol(blockw, &piv[bi*blockw], &rowbycol[bi*blockw]);
                rowbycol2colbyrow(blockw, &rowbycol[bi*blockw], &colbyrow[bi*blockw]);
            }
            if (info_)
                throw std::runtime_error("ILU failed!");
        }
        void solve(const double * const __restrict__ b, double * const __restrict__ x) const {
            // before calling solve: make sure that the 
            // block_data and sup_data pointers are still valid.
            auto y = make_unique<double[]>((this->nblocks)*(this->blockw));
            for (int bri = 0; bri < (this->nblocks); ++bri){
                for (int li = 0; li < (this->blockw); ++li){
                    double s = 0.0;
                    for (int lci = 0; lci < li; ++lci){
                        s += this->lu_get(bri, li, lci)*y[bri*(this->blockw) + lci];
                    }
                    for (int di = 1; di < (this->ndiag) + 1; ++di){
                        if (bri >= di) {
                            int ci = this->colbyrow[bri*(this->blockw) + li];
                            s += (this->sub_get(di-1, bri-di, ci) * y[(bri-di)*(this->blockw) + ci]);
                        }
                    }
                    y[bri*(this->blockw) + li] = b[bri*(this->blockw) 
                                                   + this->rowbycol[bri*(this->blockw) + li]
                                                   ] - s;
                }
            }
            for (int bri = this->nblocks - 1; bri >= 0; --bri){
                for (int li = this->blockw - 1; li >= 0; --li){
                    double s = 0.0;
                    for (int ci = li + 1; ci < (this->blockw); ++ci)
                        s += this->lu_get(bri, li, ci)*x[bri*(this->blockw) + ci];
                    for (int di = 1; di <= this->ndiag; ++di) {
                        if (bri < this->nblocks - di){
                            int ci = this->colbyrow[bri*this->blockw + li];
                            s += this->sup_get(di-1, bri, ci)*x[(bri+di)*(this->blockw) + ci];
                        }
                    }
                    x[bri*this->blockw+li] = (y[bri*(this->blockw) + li] - s)/(this->lu_get(bri, li, li));
                }
            }
        }
    };

    class BlockDiagMat {
    public:
        const int nblocks, blockw, ndiag, sub_offset, sup_offset;
        const std::size_t data_len;
        std::unique_ptr<double[]> data;
        BlockDiagMat(int nblocks, int blockw, int ndiag) :
            nblocks(nblocks), blockw(blockw), ndiag(ndiag),
            sub_offset(nblocks*blockw*blockw),
            sup_offset(nblocks*blockw*blockw + diag_store_len(nblocks, blockw, ndiag)),
            data_len(nblocks*blockw*blockw+2*diag_store_len(nblocks, blockw, ndiag)),
            data(make_unique<double[]>(data_len))
            {}
        inline void zero_out_all() noexcept {
            for (std::size_t i=0; i<(this->data_len); ++i){
                this->data[i] = 0.0;
            }
        }
        inline void zero_out_diags() noexcept {
            for (std::size_t i=(this->sub_offset); i<(this->data_len); ++i){
                this->data[i] = 0.0;
            }
        }
        inline double& block(int bi, int ri, int ci) noexcept {
            return data[bi*(this->blockw)*(this->blockw) + ci*(this->blockw) + ri];
        }
        inline double& sub(int di, int bi, int lci) noexcept {
            return data[this->sub_offset + diag_store_len(this->nblocks, this->blockw, di) + bi*(this->blockw) + lci];
        }
        inline double& sup(int di, int bi, int lci) noexcept {
            return data[this->sup_offset + diag_store_len(this->nblocks, this->blockw, di) + bi*(this->blockw) + lci];
        }

#ifdef UNIT_TEST
        double get(int ri, int ci) {
            const int bri = ri / this->blockw;
            const int bci = ci / this->blockw;
            const int lri = ri - bri*this->blockw;
            const int lci = ci - bci*this->blockw;
            if (bri == bci)
                return this->block(bri, lri, lci);
            if (std::abs(bri - bci) > ndiag)
                return 0.0;
            if (lri != lci)
                return 0.0;
            if (bri - bci > 0)
                return this->sub(bri-bci-1, bci, lci);
            return this->sup(bci-bri-1, bri, lri);
        }
#endif

        void set_to_1_minus_gamma_times_other(double gamma, BlockDiagMat &other) {
            // Scale main blocks by -gamma
            for (int bi=0; bi<this->nblocks; ++bi)
                for (int ci=0; ci<this->blockw; ++ci)
                    for (int ri=0; ri<this->blockw; ++ri)
                        this->block(bi, ri, ci) = -gamma*other.block(bi, ri, ci);

            // Add the identiy matrix
            for (int bi=0; bi<this->nblocks; ++bi)
                for (int ci=0; ci<this->blockw; ++ci)
                    this->block(bi, ci, ci) += 1;

            // Scale diagonals by -gamma
            for (int di=0; di<this->ndiag; ++di)
                for (int bi=0; bi<this->nblocks-di-1; ++bi)
                    for (int ci=0; ci<this->blockw; ++ci){
                        this->sub(di, bi, ci) = -gamma*other.sub(di, bi, ci);
                        this->sup(di, bi, ci) = -gamma*other.sup(di, bi, ci);
                    }
        }

        // The end user must assure that the underlying data is not freed.
        ILU ilu_inplace() {
            return ILU(&this->data[0],
                       &this->data[this->sub_offset],
                       &this->data[this->sup_offset],
                       this->nblocks, this->blockw, this->ndiag);
        }
        void dot_vec(const double * const vec, double * const out){
            // out need not be zeroed out before call
            const int nblocks = this->nblocks;
            const int blockw = this->blockw;
            for (int i=0; i<nblocks*blockw; ++i)
                out[i] = 0.0;
            for (int bri=0; bri<nblocks; ++bri)
                for (int lci=0; lci<blockw; ++lci)
                    for (int lri=0; lri<blockw; ++lri)
                        out[bri*blockw + lri] += vec[bri*blockw + lci]*(this->block(bri, lri, lci));
            for (int di=0; di<this->ndiag; ++di)
                for (int bi=0; bi<nblocks-di-1; ++bi)
                    for (int ci=0; ci<blockw; ++ci){
                        out[bi*blockw + ci] += this->sup(di, bi, ci)*vec[(bi+di+1)*blockw+ci];
                        out[(bi+di+1)*blockw + ci] += this->sub(di, bi, ci)*vec[bi*blockw+ci];
                    }
        }
    };

}

#if defined(WITH_BLOCK_DIAG_ILU_DGETRF)
inline int block_diag_ilu::dgetrf_square(const int dim, double * const __restrict__ a,
                                         const int lda, int * const __restrict__ ipiv){
    // Unblocked algorithm for LU decomposition of square matrices
    // employing Doolittle's algorithm with rowswaps.
    //
    // ipiv indexing starts at 1 (Fortran compability)
    if (dim == 0) return 0;

    int info = 0;
    auto A = [&](int ri, int ci) -> double& { return a[ci*lda + ri]; };
    auto swaprows = [&](int ri1, int ri2) { // this is not cache friendly
        for (int ci=0; ci<dim; ++ci){
            double temp = A(ri1, ci);
            A(ri1, ci) = A(ri2, ci);
            A(ri2, ci) = temp;
        }
    };

    for (int i=0; i<dim-1; ++i) {
        int pivrow = i;
        double absmax = std::abs(A(i, i));
        for (int j=i; j<dim; ++j) {
            // Find pivot
            double curabs = std::abs(A(j, i));
            if (curabs > absmax){
                absmax = curabs;
                pivrow = j;
            }
        }
        if ((absmax == 0) && (info == 0))
            info = pivrow+1;
        ipiv[i] = pivrow+1;
        if (pivrow != i) {
            // Swap rows
            swaprows(i, pivrow);
        }
        // Eliminate in column
        for (int ri=i+1; ri<dim; ++ri){
            A(ri, i) = A(ri, i)/A(i, i);
        }
        // Subtract from rows
        for (int ci=i+1; ci<dim; ++ci){
            for (int ri=i+1; ri<dim; ++ri){
                A(ri, ci) -= A(ri, i)*A(i, ci);
            }
        }
    }
    ipiv[dim-1] = dim;
    return info;
}
#endif

#endif
