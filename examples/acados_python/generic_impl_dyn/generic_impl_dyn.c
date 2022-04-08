#ifdef __cplusplus
extern "C" {
#endif

#include "acados/utils/external_function_generic.h"
#include "acados/utils/math.h"
#include "blasfeo_d_blas.h"
#include "blasfeo_d_aux.h"
#include "blasfeo_d_aux_ext_dep.h"
#include <math.h>

int generic_impl_dyn_fun(void **in, void **out, void *params)
{
    int ii;

    int nu = 1;
    int nx = 4;

    // extract inputs
    // 0: [x], size: nx, type: BLASFEO_DVEC_ARGS
    struct blasfeo_dvec *x = in[0];
    // 1: [xdot], size: nx, type: BLASFEO_DVEC_ARGS
    struct blasfeo_dvec_args *xdot_args = in[1];
    struct blasfeo_dvec *xdot = xdot_args->x;
    int xdoti = xdot_args->xi; // offset in vector
    // 2: [u], size: nu, type: COLMAJ
    double *u = in[2];
    // 3: [z]: BLASFEO_DVEC_ARGS, empty for this model.
    double *z = in[3];

    double theta = BLASFEO_DVECEL(x, 1);
    double v1 = BLASFEO_DVECEL(x, 2);
    double dtheta = BLASFEO_DVECEL(x, 3);

    double x1_dot = BLASFEO_DVECEL(xdot, xdoti);
    double theta_dot = BLASFEO_DVECEL(xdot, xdoti+1);
    double v1_dot = BLASFEO_DVECEL(xdot, xdoti+2);
    double dtheta_dot = BLASFEO_DVECEL(xdot, xdoti+3);

    double F = u[0];

// x
// SX([x1, theta, v1, dtheta])

// casadi f_impl
// SX(@1=-0.08, @2=(1.1-((0.1*cos(theta))*cos(theta))), [(x1_dot-v1), (theta_dot-dtheta), (v1_dot-((((((@1*sin(theta))*dtheta)*dtheta)+((0.981*cos(theta))*sin(theta)))+F)/@2)), (dtheta_dot-(((((((@1*cos(theta))*sin(theta))*dtheta)*dtheta)+(F*cos(theta)))+(10.791*sin(theta)))/(0.8*@2)))])

    // extract outputs
    // 0: [fun], size: nx1, type: BLASFEO_DVEC_ARGS
    struct blasfeo_dvec_args *f_args = out[0];
    struct blasfeo_dvec *fun = f_args->x;
    int f_off = f_args->xi; // offset in vector

    double t1 = -0.08;
    double t2 = (1.1-((0.1*cos(theta))*cos(theta)));

    double out_tmp = (x1_dot-v1);
    blasfeo_dvecse(1, out_tmp, fun, f_off+0);
    out_tmp = (theta_dot-dtheta);
    blasfeo_dvecse(1, out_tmp, fun, f_off+1);
    out_tmp = (v1_dot-((((((t1*sin(theta))*dtheta)*dtheta)+((0.981*cos(theta))*sin(theta)))+F)/t2));
    blasfeo_dvecse(1, out_tmp, fun, f_off+2);
    out_tmp = (dtheta_dot-(((((((t1*cos(theta))*sin(theta))*dtheta)*dtheta)+(F*cos(theta)))+(10.791*sin(theta)))/(0.8*t2)));
    blasfeo_dvecse(1, out_tmp, fun, f_off+3);


    return 0;
}





    // impl_dae_fun_jac_x_xdot_z = Function(fun_name, [x, xdot, u, z, p], [f_impl, jac_x, jac_xdot, jac_z])

int generic_impl_dyn_fun_jac(void **in, void **out, void *params)
{
    int ii;

    int nu = 1;
    int nx = 4;

    // extract inputs
    // 0: [x], size: nx, type: BLASFEO_DVEC_ARGS
    struct blasfeo_dvec *x = in[0];
    // 1: [xdot], size: nx, type: BLASFEO_DVEC_ARGS
    struct blasfeo_dvec_args *xdot_args = in[1];
    struct blasfeo_dvec *xdot = xdot_args->x;
    int xdoti = xdot_args->xi; // offset in vector
    // 2: [u], size: nu, type: COLMAJ
    double *u = in[2];
    // 3: [z]: BLASFEO_DVEC_ARGS, empty for this model.
    double *z = in[3];

    double theta = BLASFEO_DVECEL(x, 1);
    double v1 = BLASFEO_DVECEL(x, 2);
    double dtheta = BLASFEO_DVECEL(x, 3);

    double x1_dot = BLASFEO_DVECEL(xdot, xdoti);
    double theta_dot = BLASFEO_DVECEL(xdot, xdoti+1);
    double v1_dot = BLASFEO_DVECEL(xdot, xdoti+2);
    double dtheta_dot = BLASFEO_DVECEL(xdot, xdoti+3);

    double F = u[0];

// x
// SX([x1, theta, v1, dtheta])

// casadi f_impl
// SX(@1=-0.08, @2=(1.1-((0.1*cos(theta))*cos(theta))), [(x1_dot-v1), (theta_dot-dtheta), (v1_dot-((((((@1*sin(theta))*dtheta)*dtheta)+((0.981*cos(theta))*sin(theta)))+F)/@2)), (dtheta_dot-(((((((@1*cos(theta))*sin(theta))*dtheta)*dtheta)+(F*cos(theta)))+(10.791*sin(theta)))/(0.8*@2)))])

    // extract outputs
    // 0: [fun], size: nx1, type: BLASFEO_DVEC_ARGS
    struct blasfeo_dvec_args *f_args = out[0];
    struct blasfeo_dvec *fun = f_args->x;
    int f_off = f_args->xi; // offset in vector

    double t1 = -0.08;
    double t2 = (1.1-((0.1*cos(theta))*cos(theta)));

    double out_tmp = (x1_dot-v1);
    blasfeo_dvecse(1, out_tmp, fun, f_off+0);
    out_tmp = (theta_dot-dtheta);
    blasfeo_dvecse(1, out_tmp, fun, f_off+1);
    out_tmp = (v1_dot-((((((t1*sin(theta))*dtheta)*dtheta)+((0.981*cos(theta))*sin(theta)))+F)/t2));
    blasfeo_dvecse(1, out_tmp, fun, f_off+2);
    out_tmp = (dtheta_dot-(((((((t1*cos(theta))*sin(theta))*dtheta)*dtheta)+(F*cos(theta)))+(10.791*sin(theta)))/(0.8*t2)));
    blasfeo_dvecse(1, out_tmp, fun, f_off+3);


    // 0: jac_x
// SX(@1=-0.08, @2=cos(theta), @3=0.981, @4=cos(theta), @5=(@3*@4), @6=sin(theta), @7=sin(theta), @8=0.1, @9=(@8*@4), @10=(1.1-(@9*@4)), @11=(@1*@6), @12=(@11*dtheta), @13=((@4*(@8*@7))+(@9*@7)), @14=(@1*@4), @15=10.791, @16=0.8, @17=(@16*@10), @18=(@14*@6), @19=(@18*dtheta), @20=-1,
// [[00, 00, @20, 00],
//  [00, 00, 00, @20],
//  [00, (-((((dtheta*(dtheta*(@1*@2)))+((@5*@2)-(@6*(@3*@7))))/@10)-((((((@12*dtheta)+(@5*@6))+F)/@10)/@10)*@13))), 00, (-(((dtheta*@11)+@12)/@10))],
//  [00, (-(((((dtheta*(dtheta*((@14*@2)-(@6*(@1*@7)))))-(F*@7))+(@15*@2))/@17)-((((((@19*dtheta)+(F*@4))+(@15*@6))/@17)/@17)*(@16*@13)))), 00, (-(((dtheta*@18)+@19)/@17))]])


    struct blasfeo_dmat *jac_x = out[1];
    blasfeo_dgese(nx, nx, 0.0, jac_x, 0, 0);
    t2 = cos(theta);
    double t3 = 0.981;
    double t4 = t2;
    double t5 = t3 * t2;
    double t6 = sin(theta); // = t7
    double t8 = 0.1;
    double t9 = t8*t2;
    double t10=(1.1-(t9*t4));
    double t11=(t1*t6);
    double t12=(t11*dtheta);
    double t13=((t4*(t8*t6))+(t9*t6));
    double t14=(t1*t4);
    double t15=10.791;
    double t16=0.8;
    double t17=(t16*t10);
    double t18=(t14*t6);
    double t19=(t18*dtheta);
    double t20=-1;

    // jac_x
    blasfeo_dgein1(t20, jac_x, 0, 2);
    blasfeo_dgein1(t20, jac_x, 1, 3);

    double tmp = (-((((dtheta*(dtheta*(t1*t2)))+((t5*t2)-(t6*(t3*t6))))/t10)-((((((t12*dtheta)+(t5*t6))+F)/t10)/t10)*t13)));
    blasfeo_dgein1(tmp, jac_x, 2, 1);
    tmp = (-(((dtheta*t11)+t12)/t10));
    blasfeo_dgein1(tmp, jac_x, 2, 3);
    tmp = (-(((((dtheta*(dtheta*((t14*t2)-(t6*(t1*t6)))))-(F*t6))+(t15*t2))/t17)-((((((t19*dtheta)+(F*t4))+(t15*t6))/t17)/t17)*(t16*t13))));
    blasfeo_dgein1(tmp, jac_x, 3, 1);
    tmp = (-(((dtheta*t18)+t19)/t17));
    blasfeo_dgein1(tmp, jac_x, 3, 3);

    // Not in this fun
    // jac_u
    // SX(@1=cos(theta), @2=(1.1-((0.1*@1)*@1)), [00, 00, (-(1./@2)), (-(@1/(0.8*@2)))])

    // jac_xdot
    // SX(@1=1,
    // [[@1, 00, 00, 00],
    //  [00, @1, 00, 00],
    //  [00, 00, @1, 00],
    //  [00, 00, 00, @1]])
    struct blasfeo_dmat *jac_xdot = out[2];
    blasfeo_dgese(nx, nx, 0.0, jac_xdot, 0, 0);
    blasfeo_ddiare(nx, 1.0, jac_xdot, 0, 0);

    return 0;
}





// impl_dae_jac_x_xdot_u_z = Function(fun_name, [x, xdot, u, z, p], [jac_x, jac_xdot, jac_u, jac_z])
int generic_impl_dyn_jac(void **in, void **out, void *params)
{

    int ii;

    int nu = 1;
    int nx = 4;

    // extract inputs
    // 0: [x], size: nx, type: BLASFEO_DVEC_ARGS
    struct blasfeo_dvec *x = in[0];
    // 1: [xdot], size: nx, type: BLASFEO_DVEC_ARGS
    struct blasfeo_dvec_args *xdot_args = in[1];
    struct blasfeo_dvec *xdot = xdot_args->x;
    int xdoti = xdot_args->xi; // offset in vector
    // 2: [u], size: nu, type: COLMAJ
    double *u = in[2];
    // 3: [z]: BLASFEO_DVEC_ARGS, empty for this model.
    double *z = in[3];

    double theta = BLASFEO_DVECEL(x, 1);
    double v1 = BLASFEO_DVECEL(x, 2);
    double dtheta = BLASFEO_DVECEL(x, 3);

    double x1_dot = BLASFEO_DVECEL(xdot, xdoti);
    double theta_dot = BLASFEO_DVECEL(xdot, xdoti+1);
    double v1_dot = BLASFEO_DVECEL(xdot, xdoti+2);
    double dtheta_dot = BLASFEO_DVECEL(xdot, xdoti+3);

    double F = u[0];

// x
// SX([x1, theta, v1, dtheta])

// casadi f_impl
// SX(@1=-0.08, @2=(1.1-((0.1*cos(theta))*cos(theta))), [(x1_dot-v1), (theta_dot-dtheta), (v1_dot-((((((@1*sin(theta))*dtheta)*dtheta)+((0.981*cos(theta))*sin(theta)))+F)/@2)), (dtheta_dot-(((((((@1*cos(theta))*sin(theta))*dtheta)*dtheta)+(F*cos(theta)))+(10.791*sin(theta)))/(0.8*@2)))])

    // extract outputs
    struct blasfeo_dmat *jac_x = out[0];
    blasfeo_dgese(nx, nx, 0.0, jac_x, 0, 0);
    double t1 = -0.08;
    double t2 = cos(theta);
    double t3 = 0.981;
    double t4 = t2;
    double t5 = t3 * t2;
    double t6 = sin(theta); // = t7
    double t8 = 0.1;
    double t9 = t8*t2;
    double t10=(1.1-(t9*t4));
    double t11=(t1*t6);
    double t12=(t11*dtheta);
    double t13=((t4*(t8*t6))+(t9*t6));
    double t14=(t1*t4);
    double t15=10.791;
    double t16=0.8;
    double t17=(t16*t10);
    double t18=(t14*t6);
    double t19=(t18*dtheta);
    double t20=-1;

    // jac_x
    blasfeo_dgein1(t20, jac_x, 0, 2);
    blasfeo_dgein1(t20, jac_x, 1, 3);

    double tmp = (-((((dtheta*(dtheta*(t1*t2)))+((t5*t2)-(t6*(t3*t6))))/t10)-((((((t12*dtheta)+(t5*t6))+F)/t10)/t10)*t13)));
    blasfeo_dgein1(tmp, jac_x, 2, 1);
    tmp = (-(((dtheta*t11)+t12)/t10));
    blasfeo_dgein1(tmp, jac_x, 2, 3);
    tmp = (-(((((dtheta*(dtheta*((t14*t2)-(t6*(t1*t6)))))-(F*t6))+(t15*t2))/t17)-((((((t19*dtheta)+(F*t4))+(t15*t6))/t17)/t17)*(t16*t13))));
    blasfeo_dgein1(tmp, jac_x, 3, 1);
    tmp = (-(((dtheta*t18)+t19)/t17));
    blasfeo_dgein1(tmp, jac_x, 3, 3);


    // jac_xdot
    // SX(@1=1,
    // [[@1, 00, 00, 00],
    //  [00, @1, 00, 00],
    //  [00, 00, @1, 00],
    //  [00, 00, 00, @1]])
    struct blasfeo_dmat *jac_xdot = out[1];
    blasfeo_dgese(nx, nx, 0.0, jac_xdot, 0, 0);
    blasfeo_ddiare(nx, 1.0, jac_xdot, 0, 0);


    // jac_u
    // SX(@1=cos(theta), @2=(1.1-((0.1*@1)*@1)), [00, 00, (-(1./@2)), (-(@1/(0.8*@2)))])
    t1 = cos(theta);
    t2 = (1.1-((0.1*t1)*t1));
    struct blasfeo_dmat *jac_u = out[2];
    blasfeo_dgese(nx, nu, 0.0, jac_u, 0, 0);
    tmp = (-(1./t2));
    blasfeo_dgein1(tmp, jac_u, 2, 0);
    tmp = (-(t1/(0.8*t2)));
    blasfeo_dgein1(tmp, jac_u, 3, 0);

    return 0;
}



#ifdef __cplusplus
} /* extern "C" */
#endif
