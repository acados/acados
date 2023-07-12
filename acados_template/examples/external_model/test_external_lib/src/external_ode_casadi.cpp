#include "external_ode_casadi.h"
#ifdef __cplusplus
extern "C" {
#endif

// Declare sparsity dimensions
static bool isSparsityFilled(false);
static casadi_int x_sparsity[3] = {2, 1, 1};
static casadi_int u_sparsity[3] = {1, 1, 1};
static casadi_int xDot_sparsity[3] = {2, 1, 1};

const char* libexternal_ode_casadi_name(void){
    return "libexternal_ode_casadi";
}

int libexternal_ode_casadi(const casadi_real** arg, double** res, casadi_int*, casadi_real*, void*){

    // Return the answers
    res[0][0]   = arg[0][1];
    res[0][1]   = arg[1][0];

    return 0;
}

// IN
casadi_int libexternal_ode_casadi_n_in(void){
    return 2;
}
const char* libexternal_ode_casadi_name_in(casadi_int i){
    switch (i) {
    case 0: return "x";
    case 1: return "u";
    default: return nullptr;
    }
}
const casadi_int* libexternal_ode_casadi_sparsity_in(casadi_int i) {
    switch (i) {
    case 0: return x_sparsity;
    case 1: return u_sparsity;
    default: return nullptr;
    }
}

// OUT
casadi_int libexternal_ode_casadi_n_out(void){
    return 1;
}
const char* libexternal_ode_casadi_name_out(casadi_int i){
    switch (i) {
    case 0: return "xDot";
    default: return nullptr;
    }
}
const casadi_int* libexternal_ode_casadi_sparsity_out(casadi_int i) {
    switch (i) {
    case 0: return xDot_sparsity;
    default: return nullptr;
    }
}

int libexternal_ode_casadi_work(casadi_int *sz_arg,
                                               casadi_int* sz_res,
                                               casadi_int *sz_iw,
                                               casadi_int *sz_w) {
    if (sz_arg) *sz_arg = 3;
    if (sz_res) *sz_res = 2;
    if (sz_iw) *sz_iw = 0;
    if (sz_w) *sz_w = 0;
    return 0;
}

bool libexternal_ode_casadi_has_derivative(void)
{
    return true;
}


#ifdef __cplusplus
} /* extern "C" */
#endif
