/* The model comes from \cite{Wirsching2006} */
#include "examples/c/chain_model/chain_model.h"

extern int jac_chain_nm2(void* mem, const real_t** arg, real_t** res, int* iw, real_t* w);
extern int jac_chain_nm3(void* mem, const real_t** arg, real_t** res, int* iw, real_t* w);
extern int jac_chain_nm4(void* mem, const real_t** arg, real_t** res, int* iw, real_t* w);
extern int jac_chain_nm5(void* mem, const real_t** arg, real_t** res, int* iw, real_t* w);
extern int jac_chain_nm6(void* mem, const real_t** arg, real_t** res, int* iw, real_t* w);
extern int jac_chain_nm7(void* mem, const real_t** arg, real_t** res, int* iw, real_t* w);
extern int jac_chain_nm8(void* mem, const real_t** arg, real_t** res, int* iw, real_t* w);
extern int jac_chain_nm9(void* mem, const real_t** arg, real_t** res, int* iw, real_t* w);

// The auto-generated VDE functions from CasADi:

void VDE_fun_nm2(const real_t* in, real_t* out, 
    int (*vde)(const real_t**, real_t**, int*, real_t*, int)) {

    int_t NMF = 1;
    int_t NX = NMF*6;
    int_t NU = 3;
    const real_t* x = in;
    const real_t* Sx = in + NX;
    const real_t* Su = in + NX + NX*NX;
    const real_t* u  = in + NX + NX*(NX+NU);

    real_t* x_out = out;
    real_t* Sx_out = out + NX;
    real_t* Su_out = out + NX + NX*NX;

    int casadi_mem = 0;
    int *casadi_iw = 0;
    real_t *casadi_w = 0;

    const real_t *casadi_arg[4];
    real_t *casadi_res[3];

    casadi_arg[0] = x;
    casadi_arg[1] = Sx;
    casadi_arg[2] = Su;
    casadi_arg[3] = u;

    casadi_res[0] = x_out;
    casadi_res[1] = Sx_out;
    casadi_res[2] = Su_out;

    vde(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}

void VDE_fun_nm3(const real_t* in, real_t* out, 
    int (*vde)(const real_t**, real_t**, int*, real_t*, int)) {

    int_t NMF = 2;
    int_t NX = NMF*6;
    int_t NU = 3;
    const real_t* x = in;
    const real_t* Sx = in + NX;
    const real_t* Su = in + NX + NX*NX;
    const real_t* u  = in + NX + NX*(NX+NU);

    real_t* x_out = out;
    real_t* Sx_out = out + NX;
    real_t* Su_out = out + NX + NX*NX;

    int casadi_mem = 0;
    int *casadi_iw = 0;
    real_t *casadi_w = 0;

    const real_t *casadi_arg[4];
    real_t *casadi_res[3];

    casadi_arg[0] = x;
    casadi_arg[1] = Sx;
    casadi_arg[2] = Su;
    casadi_arg[3] = u;

    casadi_res[0] = x_out;
    casadi_res[1] = Sx_out;
    casadi_res[2] = Su_out;

    vde(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}

void VDE_fun_nm4(const real_t* in, real_t* out, 
    int (*vde)(const real_t**, real_t**, int*, real_t*, int)) {

    int_t NMF = 3;
    int_t NX = NMF*6;
    int_t NU = 3;
    const real_t* x = in;
    const real_t* Sx = in + NX;
    const real_t* Su = in + NX + NX*NX;
    const real_t* u  = in + NX + NX*(NX+NU);

    real_t* x_out = out;
    real_t* Sx_out = out + NX;
    real_t* Su_out = out + NX + NX*NX;

    int casadi_mem = 0;
    int *casadi_iw = 0;
    real_t *casadi_w = 0;

    const real_t *casadi_arg[4];
    real_t *casadi_res[3];

    casadi_arg[0] = x;
    casadi_arg[1] = Sx;
    casadi_arg[2] = Su;
    casadi_arg[3] = u;

    casadi_res[0] = x_out;
    casadi_res[1] = Sx_out;
    casadi_res[2] = Su_out;

    vde(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}

void VDE_fun_nm5(const real_t* in, real_t* out, 
    int (*vde)(const real_t**, real_t**, int*, real_t*, int)) {

    int_t NMF = 4;
    int_t NX = NMF*6;
    int_t NU = 3;
    const real_t* x = in;
    const real_t* Sx = in + NX;
    const real_t* Su = in + NX + NX*NX;
    const real_t* u  = in + NX + NX*(NX+NU);

    real_t* x_out = out;
    real_t* Sx_out = out + NX;
    real_t* Su_out = out + NX + NX*NX;

    int casadi_mem = 0;
    int *casadi_iw = 0;
    real_t *casadi_w = 0;

    const real_t *casadi_arg[4];
    real_t *casadi_res[3];

    casadi_arg[0] = x;
    casadi_arg[1] = Sx;
    casadi_arg[2] = Su;
    casadi_arg[3] = u;

    casadi_res[0] = x_out;
    casadi_res[1] = Sx_out;
    casadi_res[2] = Su_out;

    vde(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}

void VDE_fun_nm6(const real_t* in, real_t* out, 
    int (*vde)(const real_t**, real_t**, int*, real_t*, int)) {

    int_t NMF = 5;
    int_t NX = NMF*6;
    int_t NU = 3;
    const real_t* x = in;
    const real_t* Sx = in + NX;
    const real_t* Su = in + NX + NX*NX;
    const real_t* u  = in + NX + NX*(NX+NU);

    real_t* x_out = out;
    real_t* Sx_out = out + NX;
    real_t* Su_out = out + NX + NX*NX;

    int casadi_mem = 0;
    int *casadi_iw = 0;
    real_t *casadi_w = 0;

    const real_t *casadi_arg[4];
    real_t *casadi_res[3];

    casadi_arg[0] = x;
    casadi_arg[1] = Sx;
    casadi_arg[2] = Su;
    casadi_arg[3] = u;

    casadi_res[0] = x_out;
    casadi_res[1] = Sx_out;
    casadi_res[2] = Su_out;

    vde(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}

void VDE_fun_nm7(const real_t* in, real_t* out, 
    int (*vde)(const real_t**, real_t**, int*, real_t*, int)) {

    int_t NMF = 6;
    int_t NX = NMF*6;
    int_t NU = 3;
    const real_t* x = in;
    const real_t* Sx = in + NX;
    const real_t* Su = in + NX + NX*NX;
    const real_t* u  = in + NX + NX*(NX+NU);

    real_t* x_out = out;
    real_t* Sx_out = out + NX;
    real_t* Su_out = out + NX + NX*NX;

    int casadi_mem = 0;
    int *casadi_iw = 0;
    real_t *casadi_w = 0;

    const real_t *casadi_arg[4];
    real_t *casadi_res[3];

    casadi_arg[0] = x;
    casadi_arg[1] = Sx;
    casadi_arg[2] = Su;
    casadi_arg[3] = u;

    casadi_res[0] = x_out;
    casadi_res[1] = Sx_out;
    casadi_res[2] = Su_out;

    vde(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}

void VDE_fun_nm8(const real_t* in, real_t* out, 
    int (*vde)(const real_t**, real_t**, int*, real_t*, int)) {

    int_t NMF = 7;
    int_t NX = NMF*6;
    int_t NU = 3;
    const real_t* x = in;
    const real_t* Sx = in + NX;
    const real_t* Su = in + NX + NX*NX;
    const real_t* u  = in + NX + NX*(NX+NU);

    real_t* x_out = out;
    real_t* Sx_out = out + NX;
    real_t* Su_out = out + NX + NX*NX;

    int casadi_mem = 0;
    int *casadi_iw = 0;
    real_t *casadi_w = 0;

    const real_t *casadi_arg[4];
    real_t *casadi_res[3];

    casadi_arg[0] = x;
    casadi_arg[1] = Sx;
    casadi_arg[2] = Su;
    casadi_arg[3] = u;

    casadi_res[0] = x_out;
    casadi_res[1] = Sx_out;
    casadi_res[2] = Su_out;

    vde(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}

void VDE_fun_nm9(const real_t* in, real_t* out, 
    int (*vde)(const real_t**, real_t**, int*, real_t*, int)) {

    int_t NMF = 8;
    int_t NX = NMF*6;
    int_t NU = 3;
    const real_t* x = in;
    const real_t* Sx = in + NX;
    const real_t* Su = in + NX + NX*NX;
    const real_t* u  = in + NX + NX*(NX+NU);

    real_t* x_out = out;
    real_t* Sx_out = out + NX;
    real_t* Su_out = out + NX + NX*NX;

    int casadi_mem = 0;
    int *casadi_iw = 0;
    real_t *casadi_w = 0;

    const real_t *casadi_arg[4];
    real_t *casadi_res[3];

    casadi_arg[0] = x;
    casadi_arg[1] = Sx;
    casadi_arg[2] = Su;
    casadi_arg[3] = u;

    casadi_res[0] = x_out;
    casadi_res[1] = Sx_out;
    casadi_res[2] = Su_out;

    vde(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}



// The auto-generated Jacobian functions from CasADi:

void jac_fun_nm2(const real_t* in, real_t* out) {
    int_t NMF = 1;
    int_t NX = NMF*6;
    const real_t* x = in;
    const real_t* u  = in + NX;

    real_t* x_out = out;
    real_t* jac_out = out + NX;

    void *casadi_mem = 0;
    int *casadi_iw = 0;
    real_t *casadi_w = 0;

    const real_t *casadi_arg[4];
    real_t *casadi_res[3];

    casadi_arg[0] = x;
    casadi_arg[1] = u;

    casadi_res[0] = x_out;
    casadi_res[1] = jac_out;

    jac_chain_nm2(casadi_mem, casadi_arg, casadi_res, casadi_iw, casadi_w);
}

void jac_fun_nm3(const real_t* in, real_t* out) {
    int_t NMF = 2;
    int_t NX = NMF*6;
    const real_t* x = in;
    const real_t* u  = in + NX;

    real_t* x_out = out;
    real_t* jac_out = out + NX;

    void *casadi_mem = 0;
    int *casadi_iw = 0;
    real_t *casadi_w = 0;

    const real_t *casadi_arg[4];
    real_t *casadi_res[3];

    casadi_arg[0] = x;
    casadi_arg[1] = u;

    casadi_res[0] = x_out;
    casadi_res[1] = jac_out;

    jac_chain_nm3(casadi_mem, casadi_arg, casadi_res, casadi_iw, casadi_w);
}

void jac_fun_nm4(const real_t* in, real_t* out) {
    int_t NMF = 3;
    int_t NX = NMF*6;
    const real_t* x = in;
    const real_t* u  = in + NX;

    real_t* x_out = out;
    real_t* jac_out = out + NX;

    void *casadi_mem = 0;
    int *casadi_iw = 0;
    real_t *casadi_w = 0;

    const real_t *casadi_arg[4];
    real_t *casadi_res[3];

    casadi_arg[0] = x;
    casadi_arg[1] = u;

    casadi_res[0] = x_out;
    casadi_res[1] = jac_out;

    jac_chain_nm4(casadi_mem, casadi_arg, casadi_res, casadi_iw, casadi_w);
}

void jac_fun_nm5(const real_t* in, real_t* out) {
    int_t NMF = 4;
    int_t NX = NMF*6;
    const real_t* x = in;
    const real_t* u  = in + NX;

    real_t* x_out = out;
    real_t* jac_out = out + NX;

    void *casadi_mem = 0;
    int *casadi_iw = 0;
    real_t *casadi_w = 0;

    const real_t *casadi_arg[4];
    real_t *casadi_res[3];

    casadi_arg[0] = x;
    casadi_arg[1] = u;

    casadi_res[0] = x_out;
    casadi_res[1] = jac_out;

    jac_chain_nm5(casadi_mem, casadi_arg, casadi_res, casadi_iw, casadi_w);
}

void jac_fun_nm6(const real_t* in, real_t* out) {
    int_t NMF = 5;
    int_t NX = NMF*6;
    const real_t* x = in;
    const real_t* u  = in + NX;

    real_t* x_out = out;
    real_t* jac_out = out + NX;

    void *casadi_mem = 0;
    int *casadi_iw = 0;
    real_t *casadi_w = 0;

    const real_t *casadi_arg[4];
    real_t *casadi_res[3];

    casadi_arg[0] = x;
    casadi_arg[1] = u;

    casadi_res[0] = x_out;
    casadi_res[1] = jac_out;

    jac_chain_nm6(casadi_mem, casadi_arg, casadi_res, casadi_iw, casadi_w);
}

void jac_fun_nm7(const real_t* in, real_t* out) {
    int_t NMF = 6;
    int_t NX = NMF*6;
    const real_t* x = in;
    const real_t* u  = in + NX;

    real_t* x_out = out;
    real_t* jac_out = out + NX;

    void *casadi_mem = 0;
    int *casadi_iw = 0;
    real_t *casadi_w = 0;

    const real_t *casadi_arg[4];
    real_t *casadi_res[3];

    casadi_arg[0] = x;
    casadi_arg[1] = u;

    casadi_res[0] = x_out;
    casadi_res[1] = jac_out;

    jac_chain_nm7(casadi_mem, casadi_arg, casadi_res, casadi_iw, casadi_w);
}

void jac_fun_nm8(const real_t* in, real_t* out) {
    int_t NMF = 7;
    int_t NX = NMF*6;
    const real_t* x = in;
    const real_t* u  = in + NX;

    real_t* x_out = out;
    real_t* jac_out = out + NX;

    void *casadi_mem = 0;
    int *casadi_iw = 0;
    real_t *casadi_w = 0;

    const real_t *casadi_arg[4];
    real_t *casadi_res[3];

    casadi_arg[0] = x;
    casadi_arg[1] = u;

    casadi_res[0] = x_out;
    casadi_res[1] = jac_out;

    jac_chain_nm8(casadi_mem, casadi_arg, casadi_res, casadi_iw, casadi_w);
}

void jac_fun_nm9(const real_t* in, real_t* out) {
    int_t NMF = 8;
    int_t NX = NMF*6;
    const real_t* x = in;
    const real_t* u  = in + NX;

    real_t* x_out = out;
    real_t* jac_out = out + NX;

    void *casadi_mem = 0;
    int *casadi_iw = 0;
    real_t *casadi_w = 0;

    const real_t *casadi_arg[4];
    real_t *casadi_res[3];

    casadi_arg[0] = x;
    casadi_arg[1] = u;

    casadi_res[0] = x_out;
    casadi_res[1] = jac_out;

    jac_chain_nm9(casadi_mem, casadi_arg, casadi_res, casadi_iw, casadi_w);
}
