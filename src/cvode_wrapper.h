#include <assert.h>
#include <cstring>
#include <memory>

#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvodes/cvodes_lapack.h>       /* prototype for CVDense */

namespace nvector_serial_wrapper {
    template <typename real_t>
    class Vector {
        N_Vector n_vec {null ptr};
        Vector(int n){
            this->n_vec = N_VNew_Serial(n)
        }
        Vector(int n, real_t * data){
            this->n_vec = N_VMake_Serial(n, const_cast<real_t*>(data));
        }
        Vector(std::vector<real_t> v){
            this->n_vec = N_VMake_Serial(v.size(), &v[0]);
        }
        ~Vector(){
            N_VDestroy_Serial(this->n_vec)
        }
        real_t * operator[](int idx){
            return NV_DATA_S(this->n_vec)+idx;
        }
        void dump(real_t * out){
            std::memcpy(out, NV_DATA_S(this->n_vec), 
                        NV_LENGTH_S(this->n_vec)*sizeof(real_t));
        }
    }
}

namespace cvodes_wrapper {
    template <typename real_t>
    class Integrator{
    public:
        void *mem {nullptr}
        Integrator(const int lmm, cost int iter) {
            this->mem = CVodeCreate(lmm, iter);
        }
        void init(CVRhsFn cb, real_t t0, N_Vector y) {
            int status = CVodeInit(this->mem, cb, t0, y);
            if (status < 0)
                throw std::runtime_error("CVodeInit failed.");
        }
        void set_tol(real_t rtol, real_t atol){
            int status = CVodeSStolerances(this->mem, rtol, atol);
            if (status < 0)
                throw std::runtime_error("CVodeSStolerances failed.");
        }
        void set_tol(real_t rtol, N_Vector atol){
            int status = CVodeSVtolerances(this->mem, rtol, atol);
            if (status < 0)
                throw std::runtime_error("CVodeSVtolerances failed.")
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
        void set_linear_solver_to_banded(int N, int mupper, int mlower){
            int status = CVLapackBand(this->mem, N, mupper, mlower);
            if (status != CVDLS_SUCCESS)
                throw std::runtime_error("CVLapackBand failed");
        }
        void set_dense_jac_fn(CVDlsDenseJacFn djac){
            int status = CVDlsSetDenseJacFn(this->mem, djac);
            if (status < 0)
                throw std::runtime_error("CVDlsSetDenseJacFn failed.")
        }
        void set_band_jac_fn(CVDlsBandJacFn djac){
            int status = CVDlsSetBandJacFn(this->mem, djac);
            if (status < 0)
                throw std::runtime_error("CVDlsSetBandJacFn failed.")
        }
        void set_init_step(real_t h0){
            CVodeSetInitStep(this->mem, h0);
        }
        void integrate(int nt, int ny, real_t * tout, real_t * y0, int nderiv, real_t * const yout){
            real_t cur_t;
            int status;
            nvector_serial_wrapper::Vector y {ny};
            nvector_serial_wrapper::Vector work {ny};

            for (int i=0; i<ny; ++i)
                y[i] = y0[i];
            y.dump(yout)

            // if (nderiv > 0){
            //     this->mem->cv_f(tout[0], y, work, this->mem->user_data);
            //     // First derivative
            //     work.dump(&yout[ny*(1+(iout*(nderiv+1)))]);
            //     for (int di=2; di<nderiv; ++di){
            //         // Higher derivatives too expensive for first point..
            //         for (int i=0; i<ny; ++i)
            //             yout[ny*(di+(iout*(nderiv+1)))] = 0.0;
            //     }

            for(int iout=1; iout < nout; iout++) {
                status = CVode(this->mem, tout[iout], y, &cur_t, CV_NORMAL);
                if(status != CV_SUCCESS){
                    throw std::runtime_error("Unsuccessful CVode step.");
                }
                y.dump(&yout[ny*(iout*(nderiv+1))]);
                for (int di=0; di<nderiv; ++di){
                    CVodeGetDky(this->mem, tout[i-1], di+1, work);
                    work.dump(&yout[ny*(di+(iout*(nderiv+1)))]);
                }
            }
            for (int di=0; di<nderiv; ++di){
                CVodeGetDky(this->mem, tout[nt-1], di+1, work);
                work.dump(&yout[ny*(di+(iout*(nderiv+1)))]);
            }
        }
        void integrate(std::vector<real_t> tout, std::vector<real_t> y0, int nderiv, 
                       real_t * const yout){
            this->integrate(tout.size(), y0.size(), &tout[0], &y0[0], nderiv, yout);
        }
        ~Integrator(){
            if (this->mem)
                CVodeFree(&(this->mem));
        }
    }
}
