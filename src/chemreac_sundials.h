#ifndef _XFWQGDP6YJBLPH7SUYHCP7FUKQ
#define _XFWQGDP6YJBLPH7SUYHCP7FUKQ

#include <cstring> // memcpy
#include <sundials/sundials_types.h> /* def. of type realtype */
#include <nvector/nvector_serial.h>  /* N_Vector types, fcts, macros */
#include <cvode/cvode.h> /* CVODE fcts., CV_BDF, CV_ADAMS */
#include <cvode/cvode_lapack.h>
#include "chemreac.h"


namespace chemreac_sundials {

// Considerations:
// long double <=> SUNDIALS_EXTENDED_PRECISION
// double <=> SUNDIALS_DOUBLE_PRECISION
// link SUNDIALS with LAPACK.

using std::vector;
using chemreac::ReactionDiffusion;

template <typename U>
int f_cb(realtype t, N_Vector y, N_Vector ydot, void *user_data){
    U * rd = (U*)user_data;
    rd->f(t, NV_DATA_S(y), NV_DATA_S(ydot));
    return 0;
}

template <typename U>
int jac_dense_cb(long int N, realtype t, 
                 N_Vector y, N_Vector fy, DlsMat Jac, void *user_data,
                 N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
    // callback of req. signature wrapping U method.
    U * rd = (U*)user_data;
    rd->dense_jac_cmaj(t, NV_DATA_S(y), DENSE_COL(Jac, 0),
                       Jac->ldim);
    return 0;
}

template <typename U>
int jac_band_cb(long int N, long int mupper, long int mlower, realtype t, 
                N_Vector y, N_Vector fy, DlsMat Jac, void *user_data,
                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
    // callback of req. signature wrapping U method.
    U * rd = (U*)user_data;
    if (Jac->s_mu != 2*rd->n)
        throw std::runtime_error("Mismatching size of padding.");
    rd->banded_padded_jac_cmaj(t, NV_DATA_S(y), Jac->data, Jac->ldim);
    return 0;
}

template <typename T, typename U>
void direct(U * rd, 
            const vector<T> atol,
            const T rtol, const int lmm,
            const T * const __restrict__ y0, 
            const std::size_t nout,
            const T * const tout,
            T * const __restrict__ yout){
    // lmm: linear multistep method: 1: CV_ADAMS, 2: CV_BDF (cvode.h)

    int status;
    int ny = rd->n*rd->N;

    realtype cur_t; 

    // Create cvode_mem
    void *cvode_mem = nullptr;
    N_Vector interf_y = N_VMake_Serial(ny, const_cast<T*>(y0));

    cvode_mem = CVodeCreate(lmm, CV_NEWTON);
    if (cvode_mem == nullptr)
        throw std::runtime_error("CVodeCreate failed.");
    status = CVodeInit(cvode_mem, f_cb<U>, tout[0], interf_y);
    if (status < 0)
        throw std::runtime_error("CVodeInit failed.");

    // Tolerances
    if (atol.size() != 1 && atol.size() != ny)
        throw std::length_error("atol of incorrect size");
    N_Vector interf_atol = N_VNew_Serial(ny);
    for (int i=0; i < ny; ++i)
        NV_Ith_S(interf_atol, i) = 
            atol[(atol.size() == 1) ? 0 : i];
    status = CVodeSVtolerances(cvode_mem, rtol, interf_atol);
    if (status < 0)
        throw std::runtime_error("CVodeSVtolerances failed");

    status = CVodeSetUserData(cvode_mem, (void *)rd);
    if (status < 0)
        throw std::runtime_error("CVodeSetUserData failed.");

    // Setup a linear solver
    if (rd->N > 1){
        // Use a banded direct solver
        status = CVLapackBand(cvode_mem, ny, rd->n, rd->n);
        if (status != CVDLS_SUCCESS)
            throw std::runtime_error("CVLapackBand failed");
        status = CVDlsSetBandJacFn(cvode_mem, jac_band_cb<U>);
        if (status < 0)
            throw std::runtime_error("CVDlsSetBandJacFn failed.");
    } else {
        // Use a dense direct solver
        status = CVLapackDense(cvode_mem, ny);
        if (status != CVDLS_SUCCESS)
            throw std::runtime_error("CVLapackDense failed");
        status = CVDlsSetDenseJacFn(cvode_mem, jac_dense_cb<U>);
        if (status < 0)
            throw std::runtime_error("CVDlsSetDenseJacFn failed.");
    }

    memcpy(&yout[0], NV_DATA_S(interf_y), ny*sizeof(double));
    for(int iout=1; iout < nout; iout++) {
        status = CVode(cvode_mem, tout[iout], interf_y, 
                       &cur_t, CV_NORMAL);
        if(status != CV_SUCCESS)
            throw std::runtime_error("Unsuccessful CVODE step...");
        memcpy(&yout[ny*iout], NV_DATA_S(interf_y), ny*sizeof(double));
    }
    
    N_VDestroy_Serial(interf_y);
    N_VDestroy_Serial(interf_atol);
    CVodeFree(&cvode_mem);  /* Free the integrator memory */
}

} // namespace chemreac_sundials
#endif /* _XFWQGDP6YJBLPH7SUYHCP7FUKQ */
