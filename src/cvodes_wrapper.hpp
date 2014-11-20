#ifndef CHEMREAC_HRX2ZF6DAVDRVP2UH3A3BM7QLE
#define CHEMREAC_HRX2ZF6DAVDRVP2UH3A3BM7QLE

#include <assert.h>
#include <cstring>
#include <memory>

//#include <sundials/sundials_types.h> /* def. of type realtype */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvodes/cvodes.h> /* CVODE fcts., CV_BDF, CV_ADAMS */
#include <cvodes/cvodes_lapack.h>       /* prototype for CVDense */

namespace nvector_serial_wrapper {
    class Vector {
    public:
        N_Vector n_vec {nullptr};
        Vector(int n){
            this->n_vec = N_VNew_Serial(n);
        }
        Vector(int n, realtype * const data){
            this->n_vec = N_VMake_Serial(n, const_cast<realtype*>(data));
        }
        Vector(std::vector<realtype> v){
            this->n_vec = N_VMake_Serial(v.size(), &v[0]);
        }
        ~Vector(){
            N_VDestroy_Serial(this->n_vec);
        }
        realtype& operator[](int idx){
            return *(NV_DATA_S(this->n_vec)+idx);
        }
        void dump(realtype * out){
            std::memcpy(out, NV_DATA_S(this->n_vec), 
                        NV_LENGTH_S(this->n_vec)*sizeof(realtype));
        }
    };
}

namespace cvodes_wrapper {
    // Linear multistep method:
    enum class LMM : int {ADAMS=CV_ADAMS, BDF=CV_BDF};
    enum class IterType : int {NEWTON=CV_NEWTON, FUNCTIONAL=CV_FUNCTIONAL};

    class Integrator{
    public:
        void *mem {nullptr};
        Integrator(const LMM lmm, const IterType iter) {
            this->mem = CVodeCreate(static_cast<int>(lmm), static_cast<int>(iter));
        }
        void init(CVRhsFn cb, realtype t0, N_Vector y) {
            int status = CVodeInit(this->mem, cb, t0, y);
            if (status < 0)
                throw std::runtime_error("CVodeInit failed.");
        }
        void init(CVRhsFn cb, realtype t0, nvector_serial_wrapper::Vector y) {
            init(cb, t0, y.n_vec);
        }
        void init(CVRhsFn cb, realtype t0, const realtype * const y, int ny) {
            nvector_serial_wrapper::Vector yvec (ny, const_cast<realtype*>(y));
            init(cb, t0, yvec.n_vec);
        }
        void reinit(realtype t0, N_Vector y){
            CVodeReInit(this->mem, t0, y);
        }
        void reinit(realtype t0, nvector_serial_wrapper::Vector y){
            reinit(t0, y.n_vec);
        }
        void reinit(realtype t0, const realtype * const y, int ny){
            nvector_serial_wrapper::Vector yvec (ny, const_cast<realtype*>(y));
            reinit(t0, yvec.n_vec);
        }
        void set_tol(realtype rtol, realtype atol){
            int status = CVodeSStolerances(this->mem, rtol, atol);
            if (status < 0)
                throw std::runtime_error("CVodeSStolerances failed.");
        }
        void set_tol(realtype rtol, N_Vector atol){
            int status = CVodeSVtolerances(this->mem, rtol, atol);
            if (status < 0)
                throw std::runtime_error("CVodeSVtolerances failed.");
        }
        void set_tol(realtype rtol, std::vector<realtype> atol){
            nvector_serial_wrapper::Vector _atol(atol.size(), &atol[0]);
            set_tol(rtol, _atol.n_vec);
        }
        void set_user_data(void *user_data){
            int status = CVodeSetUserData(this->mem, user_data);
            if (status < 0)
                throw std::runtime_error("CVodeSetUserData failed.");
        }
        void set_linear_solver_to_dense(int ny){
            int status = CVLapackDense(this->mem, ny);
            if (status != CVDLS_SUCCESS)
                throw std::runtime_error("CVLapackDense failed");
        }
        void set_dense_jac_fn(CVDlsDenseJacFn djac){
            int status = CVDlsSetDenseJacFn(this->mem, djac);
            if (status < 0)
                throw std::runtime_error("CVDlsSetDenseJacFn failed.");
        }
        void set_linear_solver_to_banded(int N, int mupper, int mlower){
            int status = CVLapackBand(this->mem, N, mupper, mlower);
            if (status != CVDLS_SUCCESS)
                throw std::runtime_error("CVLapackBand failed");
        }
        void set_band_jac_fn(CVDlsBandJacFn djac){
            int status = CVDlsSetBandJacFn(this->mem, djac);
            if (status < 0)
                throw std::runtime_error("CVDlsSetBandJacFn failed.");
        }
        void set_init_step(realtype h0){
            CVodeSetInitStep(this->mem, h0);
        }
        void integrate(int nt, int ny, const realtype * const tout, const realtype * const y0,
                       int nderiv, realtype * const yout){
            realtype cur_t;
            int status;
            nvector_serial_wrapper::Vector y {ny};
            nvector_serial_wrapper::Vector work {ny};

            for (int i=0; i<ny; ++i)
                y[i] = y0[i];
            y.dump(yout);

            for(int iout=1; iout < nt; iout++) {
                status = CVode(this->mem, tout[iout], y.n_vec, &cur_t, CV_NORMAL);
                if(status != CV_SUCCESS){
                    throw std::runtime_error("Unsuccessful CVode step.");
                }
                y.dump(&yout[ny*(iout*(nderiv+1))]);
                for (int di=0; di<nderiv; ++di){
                    CVodeGetDky(this->mem, tout[iout-1], di+1, work.n_vec);

                    work.dump(&yout[ny*(di+(iout*(nderiv+1)))]);
                }
            }
            for (int di=0; di<nderiv; ++di){
                CVodeGetDky(this->mem, tout[nt-1], di+1, work.n_vec);
                work.dump(&yout[ny*(di+((nt-1)*(nderiv+1)))]);
            }
        }
        void integrate(const std::vector<realtype> tout, const std::vector<realtype> y0, int nderiv, 
                       realtype * const yout){
            this->integrate(tout.size(), y0.size(), &tout[0], &y0[0], nderiv, yout);
        }
        ~Integrator(){
            if (this->mem)
                CVodeFree(&(this->mem));
        }
    };

    template<class OdeSys>
    int f_cb(realtype t, N_Vector y, N_Vector ydot, void *user_data){
        OdeSys * rd = (OdeSys*)user_data;
        rd->f(t, NV_DATA_S(y), NV_DATA_S(ydot));
        return 0;
    }

    template <class OdeSys>
    int jac_dense_cb(long int N, realtype t, 
                     N_Vector y, N_Vector fy, DlsMat Jac, void *user_data,
                     N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
        // callback of req. signature wrapping OdeSys method.
        OdeSys * rd = (OdeSys*)user_data;
        rd->dense_jac_cmaj(t, NV_DATA_S(y), DENSE_COL(Jac, 0),
                           Jac->ldim);
        return 0;
    }

    template <typename OdeSys>
    int jac_band_cb(long int N, long int mupper, long int mlower, realtype t, 
                    N_Vector y, N_Vector fy, DlsMat Jac, void *user_data,
                    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
        // callback of req. signature wrapping OdeSys method.
        OdeSys * rd = (OdeSys*)user_data;
        if (Jac->s_mu != 2*rd->n)
            throw std::runtime_error("Mismatching size of padding.");
        rd->banded_padded_jac_cmaj(t, NV_DATA_S(y), Jac->data, Jac->ldim);
        return 0;
    }

    template <typename real_t, class OdeSys>
    void simple_integrate(OdeSys * rd, 
                          const std::vector<real_t> atol,
                          const real_t rtol, const int lmm,
                          const real_t * const y0, 
                          const std::size_t nout,
                          const real_t * const tout,
                          real_t * const yout){
        const int ny = rd->n*rd->N;
        Integrator integr ((lmm == CV_BDF) ? LMM::BDF : LMM::ADAMS, IterType::NEWTON);
        integr.init(f_cb<OdeSys>, tout[0], y0, ny);
        if (atol.size() == 1){
            integr.set_tol(rtol, atol[0]);
        }else{
            integr.set_tol(rtol, atol);
        }
        integr.set_user_data((void *)rd);
        if (rd->N == 1){
            integr.set_linear_solver_to_dense(rd->n);
            integr.set_dense_jac_fn(jac_dense_cb<OdeSys>);
        }else {
            integr.set_linear_solver_to_banded(ny, rd->n, rd->n);
            integr.set_band_jac_fn(jac_band_cb<OdeSys>);
        }
        integr.integrate(nout, ny, tout, y0, 0, yout);
    }
}
#endif /* CHEMREAC_HRX2ZF6DAVDRVP2UH3A3BM7QLE */
