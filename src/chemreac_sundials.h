#ifndef _XFWQGDP6YJBLPH7SUYHCP7FUKQ
#define _XFWQGDP6YJBLPH7SUYHCP7FUKQ

#include "chemreac.h"
#include <cvode/cvode_lapack.h>
#include <sundials/sundials_band.h> /* def. CVBand  */
#include <cvode/cvode.h> /* CVODE fcts., CV_BDF, CV_ADAMS */
#include <sundials/sundials_types.h> /* def. of type realtype */
#include <nvector/nvector_serial.h>  /* N_Vector types, fcts, macros */

namespace chemreac_sundials {

using chemreac::ReactionDiffusion;

int jac_cb(long int N, long int mupper, long int mlower, realtype t, 
           N_Vector y, N_Vector fy, DlsMat Jac, void * user_data,
           N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);


template <typename T>
void direct_banded(ReactionDiffusion * sys, 
                   const vector<T> atol,
                   const T rtol, const int lmm,
                   const T * const __restrict__ y0, 
                   const vector<T> tout,
                   T * const __restrict__ yout);

}; // namespace chemreac_sundials
#endif /* _XFWQGDP6YJBLPH7SUYHCP7FUKQ */
