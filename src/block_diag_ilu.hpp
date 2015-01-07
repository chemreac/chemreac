#ifndef BLOCK_DIAG_ILU_GOB3CSYR2HBHUEX4HJGA3584
#define BLOCK_DIAG_ILU_GOB3CSYR2HBHUEX4HJGA3584
#include <type_traits>
#include <utility>
#include <memory>
// C++11 source code


namespace block_diag_ilu {

    // make_unique<int[]>() only in C++14, work around:
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

    extern "C" void dgetrf_(const int* dim1, const int* dim2, double* a, int* lda, int* ipiv, int* info);

    inline void rowpiv2rowbycol(int n, const int * const piv, int * const rowbycol) {
        for (int i = 0; i < n; ++i){
            rowbycol[i] = i;
        }
        for (int i=0; i<n; ++i){
            int j = piv[i] - 1; // Fortran indexing starts with 1
            if (i != j){
                int tmp = rowbycol[j];
                rowbycol[j] = i;
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
    };

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
        int piv_get(const int idx) { return this->piv[idx]; }
        int rowbycol_get(const int idx) { return this->rowbycol[idx]; }
        int colbyrow_get(const int idx) { return this->colbyrow[idx]; }
        ILU(double * const __restrict__ block_data,
            double * const __restrict__ sub_data,
            const double * const __restrict__ sup_data,
            int nblocks, 
            int blockw, int ndiag, int ld_block_data=0)
            // block_data : column major ordering (Fortran style), 
            //     (will be overwritten)
            // sub_data : sub[0], sub[1], ... sup[ndiag-1]
            // sup_data : sup[0], sup[1], ... sup[ndiag-1]
            : block_data(block_data), sub_data(sub_data), sup_data(sup_data),
              nblocks(nblocks), blockw(blockw), ndiag(ndiag), 
              ld_block_data((ld_block_data) ? ld_block_data : blockw),
              piv(make_unique<int[]>(blockw*nblocks)),
              rowbycol(make_unique<int[]>(blockw*nblocks)),
              colbyrow(make_unique<int[]>(blockw*nblocks)) {
            int info = 0; // currently ignored
            for (int bi=0; bi<nblocks; ++bi){
                dgetrf_(&(this->blockw),
                        &(this->blockw),
                        &block_data[bi*blockw*(this->ld_block_data)],
                        &(this->ld_block_data),
                        &(this->piv[bi*blockw]),
                        &info);
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
        }
        void solve(const double * const __restrict__ b, double * const __restrict__ x) const {
            // before calling solve: make sure that the 
            // block_data and sup_data pointers are still valid.
            std::unique_ptr<double[]> y (new double[(this->nblocks)*(this->blockw)]);
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
        std::unique_ptr<double[]> data;
        BlockDiagMat(int nblocks, int blockw, int ndiag) :
            nblocks(nblocks), blockw(blockw), ndiag(ndiag),
            data(make_unique<double[]>(nblocks*blockw*blockw+2*diag_store_len(nblocks, blockw, ndiag))),
            sub_offset(nblocks*blockw*blockw), sup_offset(nblocks*blockw*blockw + diag_store_len(nblocks, blockw, ndiag)) {}
        inline double& block(int bi, int ri, int ci) {
            return data[bi*(this->blockw)*(this->blockw) + ci*(this->blockw) + ri];
        }
        inline double& sub(int di, int bi, int lci) {
            return data[this->sub_offset + diag_store_len(this->nblocks, this->blockw, di) + bi*(this->blockw) + lci];
        }
        inline double& sup(int di, int bi, int lci) {
            return data[this->sup_offset + diag_store_len(this->nblocks, this->blockw, di) + bi*(this->blockw) + lci];
        }
        void set_to_1_minus_gamma_times_other(double gamma, BlockDiagMat &other) {
            for (int bi=0; bi<this->nblocks; ++bi)
                for (int ci=0; ci<this->blockw; ++ci)
                    for (int ri=0; ri<this->blockw; ++ri)
                        this->block(bi, ri, ci) = -gamma*other.block(bi, ri, ci);

            for (int bi=0; bi<this->nblocks; ++bi)
                for (int ci=0; ci<this->blockw; ++ci)
                    this->block(bi, ci, ci) += 1;

            for (int di=0; di<this->ndiag; ++di)
                for (int bi=0; bi<this->nblocks-di-1; ++bi)
                    for (int ci=0; ci<this->blockw; ++ci){
                        this->sub(di, bi, ci) = -gamma*other.sub(di, bi, ci);
                        this->sup(di, bi, ci) = -gamma*other.sup(di, bi, ci);
                    }
        }

        // The end user must assure that the underlying data stays in scope.
        ILU ilu_inplace() {
            return ILU(&this->data[0],
                       &this->data[this->sub_offset],
                       &this->data[this->sup_offset],
                       this->nblocks, this->blockw, this->ndiag);
        }
        void dot_vec(const double * const vec, double * const out){
            const int nblocks = this->nblocks;
            const int blockw = this->blockw;
            for (int i=0; i<nblocks*blockw; ++i)
                out[i] = 0.0;
            for (int bi=0; bi<nblocks; ++bi)
                for (int ci=0; ci<blockw; ++ci)
                    for (int ri=0; ri<blockw; ++ri)
                        out[bi*blockw + ri] += vec[bi*blockw + ci]*(this->block(bi, ri, ci));
            for (int di=0; di<this->ndiag; ++di)
                for (int bi=0; bi<nblocks-di-1; ++bi)
                    for (int ci=0; ci<blockw; ++ci){
                        out[bi*blockw + ci] += this->sup(di, bi, ci)*vec[(bi+di+1)*blockw+ci];
                        out[(bi+di+1)*blockw + ci] += this->sub(di, bi, ci)*vec[bi*blockw+ci];
                    }
        }
    };

};
#endif
