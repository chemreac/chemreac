#include "chemreac_sundials.h"

// Considerations:
// long double <=> SUNDIALS_EXTENDED_PRECISION
// double <=> SUNDIALS_DOUBLE_PRECISION
// link SUNDIALS with LAPACK.

namespace chemreac_sundials {

int jac_cb(long int N, long int mupper, long int mlower, realtype t, 
           N_Vector y, N_Vector fy, DlsMat Jac, void * user_data,
           N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
    // callback of req. signature wrapping ReactionDiffusion method.
    ReactionDiffusion * sys = (ReactionDiffusion*)user_data;
    if (Jac->s_mu != sys->n)
        throw "Mismatching size of padding.";
    sys->banded_padded_jac_cmaj(t, NV_DATA_S(y), BAND_COL(Jac, 0),
                                Jac->ldim);
    return 0;
}

template <typename T>
void direct_banded(ReactionDiffusion * sys,
                   const vector<T> atol, 
                   const T rtol, const int lmm,
                   const T * const __restrict__ y0, 
                   const vector<T> tout,
                   T * const __restrict__ yout)
{
    // lmm: linear multistep method: 1: CV_ADAMS, 2: CV_BDF (cvode.h)

    int status;
    int ny = sys->n*sys->N;

    realtype cur_t; 

    // Create cvode_mem
    void *cvode_mem = nullptr;
    N_Vector interf_y = N_VMake_Serial(ny, y0);
    cvode_mem = CVodeCreate(lmm, CV_NEWTON);
    if (cvode_mem == nullptr)
        throw std::bad_alloc("CVodeCreate failed.");
    status = CVodeInit(cvode_mem, sys->f, tout[0], y0);
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

    status = CVodeSetUserData(cvode_mem, (void *)&sys);
    if (status < 0)
        throw std::runtime_error("CVodeSetUserData failed.");

    status = CVLapackBand(cvode_mem, ny, sys->n, sys->n);
    if (status < 0)
        throw std::runtime_error("CVLapackBand failed");

    status = CVDlsSetBandJacFn(cvode_mem, jac_cb);
    if (status < 0)
        throw std::runtime_error("CVDlsSetBandJacFn failed.");

    memcpy(&yout[0], NV_DATA_S(interf_y), ny);
    for(int iout=1; iout < tout.size(); iout++) {
        status = CVode(cvode_mem, tout[iout], interf_y, 
                       &cur_t, CV_NORMAL);
        if(status != CV_SUCCESS)
            throw std::runtime_error("Unsuccessful CVODE step...");
        memcpy(&yout[ny*iout], NV_DATA_S(interf_y), ny);
    }
    
    N_VDestroy_Serial(interf_y);
    N_VDestroy_Serial(interf_atol);
    CVodeFree(&cvode_mem);  /* Free the integrator memory */
}

}; // namespace chemreac_sundials
