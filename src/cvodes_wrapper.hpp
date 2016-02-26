#ifndef CHEMREAC_HRX2ZF6DAVDRVP2UH3A3BM7QLE
#define CHEMREAC_HRX2ZF6DAVDRVP2UH3A3BM7QLE

#include <assert.h>
#include <cstring>
#include <memory>

//#include <sundials/sundials_types.h> /* def. of type realtype */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvodes/cvodes_spils.h>
#include <cvodes/cvodes_spgmr.h>
#include <cvodes/cvodes_spbcgs.h>
#include <cvodes/cvodes_sptfqmr.h>
#include <cvodes/cvodes.h> /* CVODE fcts., CV_BDF, CV_ADAMS */
#include <cvodes/cvodes_lapack.h>       /* prototype for CVDense */

#include <iostream> // DEBUG

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
    // Wrapper for Revision 1.34 of cvodes

    // Linear multistep method:
    enum class LMM : int {ADAMS=CV_ADAMS, BDF=CV_BDF};
    enum class IterType : int {NEWTON=CV_NEWTON, FUNCTIONAL=CV_FUNCTIONAL};
    enum class IterLinSolEnum : int {GMRES=1, BICGSTAB=2, TFQMR=3};
    enum class PrecType : int {NONE=PREC_NONE, LEFT=PREC_LEFT,
            RIGHT=PREC_RIGHT, BOTH=PREC_BOTH};

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
        void init(CVRhsFn cb, realtype t0, nvector_serial_wrapper::Vector &y) {
            init(cb, t0, y.n_vec);
        }
        void init(CVRhsFn cb, realtype t0, const realtype * const y, int ny) {
            nvector_serial_wrapper::Vector yvec (ny, const_cast<realtype*>(y));
            init(cb, t0, yvec.n_vec);
            // it is ok that yvec is destructed here
            // (see CVodeInit in cvodes.c which at L843 calls cvAllocVectors (L3790) which )
        }
        void reinit(realtype t0, N_Vector y){
            CVodeReInit(this->mem, t0, y);
        }
        void reinit(realtype t0, nvector_serial_wrapper::Vector &y){
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
        // Iterative Linear solvers
        void set_linear_solver_to_iterative(IterLinSolEnum solver, int maxl=0){
            int flag;
            switch (solver) {
            case IterLinSolEnum::GMRES:
                flag = CVSpgmr(this->mem, (int)PrecType::LEFT, maxl);
                CVSpilsSetGSType(this->mem, MODIFIED_GS); // FIX THIS
                break;
            case IterLinSolEnum::BICGSTAB:
                flag = CVSpbcg(this->mem, (int)PrecType::LEFT, maxl);
                break;
            case IterLinSolEnum::TFQMR:
                flag = CVSptfqmr(this->mem, (int)PrecType::LEFT, maxl);
                break;
            }
            switch (flag){
            case CVSPILS_SUCCESS:
                break;
            case CVSPILS_MEM_NULL:
                throw std::runtime_error("set_linear_solver_to_iterative failed (cvode_mem is NULL)");
            case CVSPILS_ILL_INPUT:
                throw std::runtime_error("PREC_LEFT invalid.");
            case CVSPILS_MEM_FAIL:
                throw std::runtime_error("Memory allocation request failed.");
            }
        }
        void cvspils_check_flag(int flag, bool check_ill_input=false) {
            switch (flag){
            case CVSPILS_SUCCESS:
                break;
            case CVSPILS_MEM_NULL:
                throw std::runtime_error("cvode_mem is NULL");
            case CVSPILS_LMEM_NULL:
                throw std::runtime_error("CVSPILS linear solver has not been initialized)");
            }
            if ((check_ill_input) && (flag == CVSPILS_ILL_INPUT))
                throw std::runtime_error("Bad input.");
        }
        void set_jac_times_vec_fn(CVSpilsJacTimesVecFn jac_times_vec_fn){
            int flag = CVSpilsSetJacTimesVecFn(this->mem, jac_times_vec_fn);
            this->cvspils_check_flag(flag);
        }
        void set_preconditioner(CVSpilsPrecSetupFn setup_fn, CVSpilsPrecSolveFn solve_fn){
            int flag = CVSpilsSetPreconditioner(this->mem, setup_fn, solve_fn);
            this->cvspils_check_flag(flag);
        }
        void set_iter_eps_lin(realtype delta){
            int flag = CVSpilsSetEpsLin(this->mem, delta);
            this->cvspils_check_flag(flag, true);
        }
        void set_prec_type(PrecType pretyp){
            int flag = CVSpilsSetPrecType(this->mem, (int)pretyp);
            this->cvspils_check_flag(flag, true);
        }

        void set_init_step(realtype h0){
            CVodeSetInitStep(this->mem, h0);
        }

        long int get_n_lin_iters(){
            long int res=0;
            int flag;
            flag = CVSpilsGetNumLinIters(this->mem, &res);
            this->cvspils_check_flag(flag);
            return res;
        }

        long int get_n_prec_evals(){
            long int res=0;
            int flag;
            flag = CVSpilsGetNumPrecEvals(this->mem, &res);
            this->cvspils_check_flag(flag);
            return res;
        }

        long int get_n_prec_solves(){
            long int res=0;
            int flag;
            flag = CVSpilsGetNumPrecSolves(this->mem, &res);
            this->cvspils_check_flag(flag);
            return res;
        }

        long int get_n_conv_fails(){
            long int res=0;
            int flag;
            flag = CVSpilsGetNumConvFails(this->mem, &res);
            this->cvspils_check_flag(flag);
            return res;
        }

        long int get_n_jac_times_evals(){
            long int res=0;
            int flag;
            flag = CVSpilsGetNumJtimesEvals(this->mem, &res);
            this->cvspils_check_flag(flag);
            return res;
        }

        long int get_n_iter_rhs(){
            long int res=0;
            int flag;
            flag = CVSpilsGetNumRhsEvals(this->mem, &res);
            this->cvspils_check_flag(flag);
            return res;
        }


        void cv_check_flag(int flag) {
            switch (flag){
            case CV_SUCCESS:
                break;
            case CV_MEM_NULL:
                throw std::runtime_error("cvode_mem is NULL");
            }
        }

        long int get_n_steps(){
            long int res=0;
            int flag = CVodeGetNumSteps(this->mem, &res);
            cv_check_flag(flag);
            return res;
        }

        long int get_n_rhs_evals(){
            long int res=0;
            int flag = CVodeGetNumRhsEvals(this->mem, &res);
            cv_check_flag(flag);
            return res;
        }

        long int get_n_lin_solv_setups(){
            long int res=0;
            int flag = CVodeGetNumLinSolvSetups(this->mem, &res);
            cv_check_flag(flag);
            return res;
        }

        long int get_n_err_test_fails(){
            long int res=0;
            int flag = CVodeGetNumErrTestFails(this->mem, &res);
            cv_check_flag(flag);
            return res;
        }

        long int get_n_nonlin_solv_iters(){
            long int res=0;
            int flag = CVodeGetNumNonlinSolvIters(this->mem, &res);
            cv_check_flag(flag);
            return res;
        }

        long int get_n_nonlin_solv_conv_fails(){
            long int res=0;
            int flag = CVodeGetNumNonlinSolvConvFails(this->mem, &res);
            cv_check_flag(flag);
            return res;
        }

        void cvdls_check_flag(int flag) {
            switch (flag){
            case CVDLS_SUCCESS:
                break;
            case CVDLS_MEM_NULL:
                throw std::runtime_error("cvode_mem is NULL");
            case CVDLS_LMEM_NULL:
                throw std::runtime_error("CVDLS linear solver has not been initialized)");
            }
        }
        long int get_n_dls_jac_evals(){
            long int res=0;
            int flag = CVDlsGetNumJacEvals(this->mem, &res);
            cvdls_check_flag(flag);
            return res;
        }

        long int get_n_dls_rhs_evals(){
            long int res=0;
            int flag = CVDlsGetNumRhsEvals(this->mem, &res);
            cvdls_check_flag(flag);
            return res;
        }

        void get_dky(realtype t, int k, nvector_serial_wrapper::Vector &dky) {
            int flag = CVodeGetDky(this->mem, t, k, dky.n_vec);
            switch(flag){
            case CV_SUCCESS:
                // CVodeGetDky succeeded.
                break;
            case CV_BAD_K:
                throw std::runtime_error("CVodeGetDky failed with (invalid k)");
                break;
            case CV_BAD_T:
                throw std::runtime_error("CVodeGetDky failed with (invalid t)");
                break;
            case CV_BAD_DKY:
                throw std::runtime_error("CVodeGetDky failed with (dky.n_vec was NULL)");
                break;
            case CV_MEM_NULL:
                throw std::runtime_error("CVodeGetDky failed with (cvode_mem was NULL)");
                break;
            }
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

            bool early_exit = false;
            for(int iout=1; (iout < nt) && (!early_exit); iout++) {
                status = CVode(this->mem, tout[iout], y.n_vec, &cur_t, CV_NORMAL);
                if(status != CV_SUCCESS){
                    early_exit = true;
                    //throw std::runtime_error("Unsuccessful CVodes step.");
                }
                y.dump(&yout[ny*(iout*(nderiv+1))]);
                for (int di=0; di<nderiv; ++di){
                    this->get_dky(tout[iout-1], di+1, work);
                    work.dump(&yout[ny*(di+(iout*(nderiv+1)))]);
                }
            }
            for (int di=0; di<nderiv; ++di){
                this->get_dky(tout[nt-1], di+1, work);
                work.dump(&yout[ny*(di+((nt-1)*(nderiv+1)))]);
            }
        }

        void integrate(const std::vector<realtype> tout, const std::vector<realtype> y0,
                       int nderiv, realtype * const yout){
            this->integrate(tout.size(), y0.size(), &tout[0], &y0[0], nderiv, yout);
        }

        ~Integrator(){
            if (this->mem)
                CVodeFree(&(this->mem));
        }
    };


    template<class OdeSys>
    int f_cb(realtype t, N_Vector y, N_Vector ydot, void *user_data){
        OdeSys * odesys = (OdeSys*)user_data;
        odesys->f(t, NV_DATA_S(y), NV_DATA_S(ydot));
        return 0;
    }

    template<class T> void ignore( const T& ) { } // ignore compiler warnings about unused parameter

    template <class OdeSys>
    int jac_dense_cb(long int N, realtype t,
                     N_Vector y, N_Vector fy, DlsMat Jac, void *user_data,
                     N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
        // callback of req. signature wrapping OdeSys method.
        ignore(N); ignore(tmp1); ignore(tmp2); ignore(tmp3);
        OdeSys * odesys = (OdeSys*)user_data;
        odesys->dense_jac_cmaj(t, NV_DATA_S(y), NV_DATA_S(fy), DENSE_COL(Jac, 0),
                           Jac->ldim);
        return 0;
    }

    template <typename OdeSys>
    int jac_band_cb(long int N, long int mupper, long int mlower, realtype t,
                    N_Vector y, N_Vector fy, DlsMat Jac, void *user_data,
                    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
        // callback of req. signature wrapping OdeSys method.
        ignore(N); ignore(mupper); ignore(mlower); ignore(tmp1); ignore(tmp2); ignore(tmp3);
        OdeSys * odesys = (OdeSys*)user_data;
        if (Jac->s_mu != 2*(odesys->n))
            throw std::runtime_error("Mismatching size of padding.");
        odesys->banded_padded_jac_cmaj(t, NV_DATA_S(y), NV_DATA_S(fy), Jac->data, Jac->ldim);
        return 0;
    }


    template <typename OdeSys>
    int jac_times_vec_cb(N_Vector v, N_Vector Jv, realtype t, N_Vector y,
                         N_Vector fy, void *user_data, N_Vector tmp){
        // callback of req. signature wrapping OdeSys method.
        ignore(tmp);
        OdeSys * odesys = (OdeSys*)user_data;
        odesys->jac_times_vec(NV_DATA_S(v), NV_DATA_S(Jv), t, NV_DATA_S(y), NV_DATA_S(fy));
        return 0;
    }

    template <typename OdeSys>
    int jac_prec_solve_cb(realtype t, N_Vector y, N_Vector fy, N_Vector r,
                          N_Vector z, realtype gamma, realtype delta, int lr,
                          void *user_data, N_Vector tmp){
        // callback of req. signature wrapping OdeSys method.
        //std::cout << "in jac_prec_solve_cb() with lr=" << lr << std::endl; // DEBUG
        ignore(tmp); ignore(delta);
        OdeSys * odesys = (OdeSys*)user_data;
        if (lr != 1)
            throw std::runtime_error("Only left preconditioning implemented.");
        odesys->prec_solve_left(t, NV_DATA_S(y), NV_DATA_S(fy), NV_DATA_S(r),
                                NV_DATA_S(z), gamma);
        return 0; // Direct solver give no hint on success, hence report success.
    }

    template <typename OdeSys>
    int prec_setup_cb(realtype t, N_Vector y, N_Vector fy, booleantype jok,
                      booleantype *jcurPtr, realtype gamma, void *user_data,
                      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
        // callback of req. signature wrapping OdeSys method.
        //std::cout << "in prec_setup_cb()" << std::endl; // DEBUG
        ignore(tmp1); ignore(tmp2); ignore(tmp3);
        OdeSys * odesys = (OdeSys*)user_data;
        bool jac_recomputed;
        odesys->prec_setup(t, NV_DATA_S(y), NV_DATA_S(fy),
                           (jok == TRUE) ? true : false,
                           jac_recomputed, gamma);
        (*jcurPtr) = (jac_recomputed) ? TRUE : FALSE;
        return 0;
    }

    template <typename real_t, class OdeSys>
    void simple_integrate(OdeSys * rd,
                          const std::vector<real_t> atol,
                          const real_t rtol, const int lmm,
                          const real_t * const y0,
                          const std::size_t nout,
                          const real_t * const tout,
                          real_t * const yout,
                          bool with_jacobian=true,
                          int iterative=0){
        // iterative == 0 => direct (Newton)
        // iterative == 1 => iterative (GMRES)
        // iterative == 2 => iterative (BiCGStab)
        // iterative == 3 => iterative (TFQMR)
        const int ny = rd->n*rd->N;
        Integrator integr {(lmm == CV_BDF) ? LMM::BDF : LMM::ADAMS,
                IterType::NEWTON};
                //(iterative) ? IterType::FUNCTIONAL : IterType::NEWTON};
        integr.set_user_data((void *)rd);
        integr.init(f_cb<OdeSys>, tout[0], y0, ny);
        if (atol.size() == 1){
            integr.set_tol(rtol, atol[0]);
        }else{
            integr.set_tol(rtol, atol);
        }
        if (rd->N == 1){
            if (iterative)
                throw std::runtime_error("Iterative solution not implemented for N==1");
            integr.set_linear_solver_to_dense(rd->n);
            if (with_jacobian)
                integr.set_dense_jac_fn(jac_dense_cb<OdeSys>);
        }else {
            if (iterative){
                switch (iterative) {
                case 1:
                    integr.set_linear_solver_to_iterative(IterLinSolEnum::GMRES); break;
                case 2:
                    integr.set_linear_solver_to_iterative(IterLinSolEnum::BICGSTAB); break;
                case 3:
                    integr.set_linear_solver_to_iterative(IterLinSolEnum::TFQMR); break;
                }
                integr.set_prec_type(PrecType::LEFT);
                integr.set_iter_eps_lin(0); // 0 => default.
                integr.set_jac_times_vec_fn(jac_times_vec_cb<OdeSys>);
                integr.set_preconditioner(prec_setup_cb<OdeSys>,
                                          jac_prec_solve_cb<OdeSys>);
                // integr.set_gram_schmidt_type() // GMRES
                // integr.set_krylov_max_len()  // BiCGStab, TFQMR
            } else {
                integr.set_linear_solver_to_banded(ny, rd->n, rd->n);
                if (with_jacobian)
                    integr.set_band_jac_fn(jac_band_cb<OdeSys>);
            }
        }
        integr.integrate(nout, ny, tout, y0, 0, yout);
        // BEGIN DEBUG
        std::cout << "n_steps=" << integr.get_n_steps() << std::endl;
        std::cout << "n_rhs_evals=" << integr.get_n_rhs_evals() << std::endl;
        std::cout << "n_lin_solv_setups=" << integr.get_n_lin_solv_setups() << std::endl;
        std::cout << "n_err_test_fails=" << integr.get_n_err_test_fails() << std::endl;
        std::cout << "n_nonlin_solv_iters=" << integr.get_n_nonlin_solv_iters() << std::endl;
        std::cout << "n_nonlin_solv_conv_fails=" << integr.get_n_nonlin_solv_conv_fails() << std::endl;
        if (iterative) {
            std::cout << "Krylov specific:" << std::endl;
            std::cout << "  n_lin_iters=" << integr.get_n_lin_iters() << std::endl;
            std::cout << "  n_prec_evals=" << integr.get_n_prec_evals() << std::endl;
            std::cout << "  n_prec_solves=" << integr.get_n_prec_solves() << std::endl;
            std::cout << "  n_conv_fails=" << integr.get_n_conv_fails() << std::endl;
            std::cout << "  n_jac_times_evals=" << integr.get_n_jac_times_evals() << std::endl;
            std::cout << "  n_iter_rhs=" << integr.get_n_iter_rhs() << std::endl;
            std::cout.flush();
        } else {
            std::cout << "Dense linear solver specific:" << std::endl;
            std::cout << "  n_dls_jac_evals=" << integr.get_n_dls_jac_evals() << std::endl;
            std::cout << "  n_dls_rhs_evals=" << integr.get_n_dls_rhs_evals() << std::endl;
        }
        // END DEBUG
    }
}
#endif /* CHEMREAC_HRX2ZF6DAVDRVP2UH3A3BM7QLE */
