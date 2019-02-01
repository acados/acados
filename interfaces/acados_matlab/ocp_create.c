// system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
// acados
#include "acados_c/ocp_nlp_interface.h"
// mex
#include "mex.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
	{

//	mexPrintf("\nin sim_create\n");

	// sizeof(long long) == sizeof(void *) = 64 !!!
	long long *l_ptr;
	int *i_ptr;
	int *tmp_idx;
	double *d_ptr;

	int ii, jj, idx;



	/* RHS */

	// model

	// dims
	double T;		bool set_T = false;
	int nx;			bool set_nx = false;
	int nu;			bool set_nu = false;
	int ny;			bool set_ny = false;
	int ny_e;		bool set_ny_e = false;
	int nbx = 0;	bool set_nbx = false;
	int nbu;		bool set_nbu = false;
	int ng;			bool set_ng = false;
	int ng_e;		bool set_ng_e = false;
	int nh;			bool set_nh = false;
	int nh_e;		bool set_nh_e = false;
	int ns = 0;		bool set_ns = false;
	int ns_e = 0;	bool set_ns_e = false;
	int nsbu = 0;	bool set_nsbu = false;
	int nsbx = 0;	bool set_nsbx = false;
	int nsg = 0;	bool set_nsg = false;
	int nsg_e = 0;	bool set_nsg_e = false;
	int nsh = 0;	bool set_nsh = false;
	int nsh_e = 0;	bool set_nsh_e = false;
	// cost
	char *cost_type;
	char *cost_e_type;
	double *Vu;		bool set_Vu = false;
	double *Vx;		bool set_Vx = false;
	double *Vx_e;	bool set_Vx_e = false;
	double *W;		bool set_W = false;
	double *W_e;	bool set_W_e = false;
	double *yr;		bool set_yr = false;
	double *yr_e;	bool set_yr_e = false;
	double *Z;		bool set_Z = false;
	double *Z_e;	bool set_Z_e = false;
	double *Zl;		bool set_Zl = false;
	double *Zl_e;	bool set_Zl_e = false;
	double *Zu;		bool set_Zu = false;
	double *Zu_e;	bool set_Zu_e = false;
	double *z;		bool set_z = false;
	double *z_e;	bool set_z_e = false;
	double *zl;		bool set_zl = false;
	double *zl_e;	bool set_zl_e = false;
	double *zu;		bool set_zu = false;
	double *zu_e;	bool set_zu_e = false;
	// constraints
	char *constr_type;
	double *x0;		bool set_x0 = false;
	double *Jbx;	bool set_Jbx = false;
	double *lbx;	bool set_lbx = false;
	double *ubx;	bool set_ubx = false;
	double *Jbu;	bool set_Jbu = false;
	double *lbu;	bool set_lbu = false;
	double *ubu;	bool set_ubu = false;
	double *C;		bool set_C = false;
	double *D;		bool set_D = false;
	double *lg;		bool set_lg = false;
	double *ug;		bool set_ug = false;
	double *C_e;	bool set_C_e = false;
	double *lg_e;	bool set_lg_e = false;
	double *ug_e;	bool set_ug_e = false;
	double *lh;		bool set_lh = false;
	double *uh;		bool set_uh = false;
	double *lh_e;	bool set_lh_e = false;
	double *uh_e;	bool set_uh_e = false;
	double *Jsbu;	bool set_Jsbu = false;
//	double *lsbu;	bool set_lsbu = false;
//	double *usbu;	bool set_usbu = false;
	double *Jsbx;	bool set_Jsbx = false;
//	double *lsbx;	bool set_lsbx = false;
//	double *usbx;	bool set_usbx = false;
	double *Jsg;	bool set_Jsg = false;
//	double *lsg;	bool set_lsg = false;
//	double *usg;	bool set_usg = false;
	double *Jsg_e;	bool set_Jsg_e = false;
//	double *lsg_e;	bool set_lsg_e = false;
//	double *usg_e;	bool set_usg_e = false;
	double *Jsh;	bool set_Jsh = false;
//	double *lsh;	bool set_lsh = false;
//	double *ush;	bool set_ush = false;
	double *Jsh_e;	bool set_Jsh_e = false;
//	double *lsh_e;	bool set_lsh_e = false;
//	double *ush_e;	bool set_ush_e = false;
	// dynamics
	char *dyn_type;
	// trajectory initialization
	double *x_init; bool set_x_init = false;
	double *u_init; bool set_u_init = false;

	// dims
	// T
	if(mxGetField( prhs[0], 0, "T" )!=NULL)
		{
		set_T = true;
		T = mxGetScalar( mxGetField( prhs[0], 0, "T" ) );
		}
	else
		{
		mexPrintf("\nerror: ocp_create: T not set!\n");
		return;
		}
	// nx
	if(mxGetField( prhs[0], 0, "nx" )!=NULL)
		{
		set_nx = true;
		nx = mxGetScalar( mxGetField( prhs[0], 0, "nx" ) );
		}
	else
		{
		mexPrintf("\nerror: ocp_create: nx not set!\n");
		return;
		}
	// nu
	if(mxGetField( prhs[0], 0, "nu" )!=NULL)
		{
		set_nu = true;
		nu = mxGetScalar( mxGetField( prhs[0], 0, "nu" ) );
		}
	else
		{
		mexPrintf("\nerror: ocp_create: nu not set!\n");
		return;
		}
	// ny
	if(mxGetField( prhs[0], 0, "ny" )!=NULL)
		{
		set_ny = true;
		ny = mxGetScalar( mxGetField( prhs[0], 0, "ny" ) );
		}
	// ny_e
	if(mxGetField( prhs[0], 0, "ny_e" )!=NULL)
		{
		set_ny_e = true;
		ny_e = mxGetScalar( mxGetField( prhs[0], 0, "ny_e" ) );
		}
	// nbx
	if(mxGetField( prhs[0], 0, "nbx" )!=NULL)
		{
		set_nbx = true;
		nbx = mxGetScalar( mxGetField( prhs[0], 0, "nbx" ) );
		}
	// nbu
	if(mxGetField( prhs[0], 0, "nbu" )!=NULL)
		{
		set_nbu = true;
		nbu = mxGetScalar( mxGetField( prhs[0], 0, "nbu" ) );
		}
	// ng
	if(mxGetField( prhs[0], 0, "ng" )!=NULL)
		{
		set_ng = true;
		ng = mxGetScalar( mxGetField( prhs[0], 0, "ng" ) );
		}
	// ng_e
	if(mxGetField( prhs[0], 0, "ng_e" )!=NULL)
		{
		set_ng_e = true;
		ng_e = mxGetScalar( mxGetField( prhs[0], 0, "ng_e" ) );
		}
	// nh
	if(mxGetField( prhs[0], 0, "nh" )!=NULL)
		{
		set_nh = true;
		nh = mxGetScalar( mxGetField( prhs[0], 0, "nh" ) );
		}
	// nh_e
	if(mxGetField( prhs[0], 0, "nh_e" )!=NULL)
		{
		set_nh_e = true;
		nh_e = mxGetScalar( mxGetField( prhs[0], 0, "nh_e" ) );
		}
	// ns
	if(mxGetField( prhs[0], 0, "ns" )!=NULL)
		{
		set_ns = true;
		ns = mxGetScalar( mxGetField( prhs[0], 0, "ns" ) );
		}
	// ns_e
	if(mxGetField( prhs[0], 0, "ns_e" )!=NULL)
		{
		set_ns_e = true;
		ns_e = mxGetScalar( mxGetField( prhs[0], 0, "ns_e" ) );
		}
	// nsbu
	if(mxGetField( prhs[0], 0, "nsbu" )!=NULL)
		{
		set_nsbu = true;
		nsbu = mxGetScalar( mxGetField( prhs[0], 0, "nsbu" ) );
		}
	// nsbx
	if(mxGetField( prhs[0], 0, "nsbx" )!=NULL)
		{
		set_nsbx = true;
		nsbx = mxGetScalar( mxGetField( prhs[0], 0, "nsbx" ) );
		}
	// nsg
	if(mxGetField( prhs[0], 0, "nsg" )!=NULL)
		{
		set_nsg = true;
		nsg = mxGetScalar( mxGetField( prhs[0], 0, "nsg" ) );
		}
	// nsg_e
	if(mxGetField( prhs[0], 0, "nsg_e" )!=NULL)
		{
		set_nsg_e = true;
		nsg_e = mxGetScalar( mxGetField( prhs[0], 0, "nsg_e" ) );
		}
	// nsh
	if(mxGetField( prhs[0], 0, "nsh" )!=NULL)
		{
		set_nsh = true;
		nsh = mxGetScalar( mxGetField( prhs[0], 0, "nsh" ) );
		}
	// nsh_e
	if(mxGetField( prhs[0], 0, "nsh_e" )!=NULL)
		{
		set_nsh_e = true;
		nsh_e = mxGetScalar( mxGetField( prhs[0], 0, "nsh_e" ) );
		}
	// cost
	// cost type
	if(mxGetField( prhs[0], 0, "cost_type" )!=NULL)
		{
		cost_type = mxArrayToString( mxGetField( prhs[0], 0, "cost_type" ) );
		}
	// cost e_type
	if(mxGetField( prhs[0], 0, "cost_e_type" )!=NULL)
		{
		cost_e_type = mxArrayToString( mxGetField( prhs[0], 0, "cost_e_type" ) );
		}
	// Vu
	if(mxGetField( prhs[0], 0, "Vu" )!=NULL)
		{
		set_Vu = true;
		Vu = mxGetPr( mxGetField( prhs[0], 0, "Vu" ) );
		}
	// Vx
	if(mxGetField( prhs[0], 0, "Vx" )!=NULL)
		{
		set_Vx = true;
		Vx = mxGetPr( mxGetField( prhs[0], 0, "Vx" ) );
		}
	// Vx_e
	if(mxGetField( prhs[0], 0, "Vx_e" )!=NULL)
		{
		set_Vx_e = true;
		Vx_e = mxGetPr( mxGetField( prhs[0], 0, "Vx_e" ) );
		}
	// W
	if(mxGetField( prhs[0], 0, "W" )!=NULL)
		{
		set_W = true;
		W = mxGetPr( mxGetField( prhs[0], 0, "W" ) );
		}
	// W_e
	if(mxGetField( prhs[0], 0, "W_e" )!=NULL)
		{
		set_W_e = true;
		W_e = mxGetPr( mxGetField( prhs[0], 0, "W_e" ) );
		}
	// yr
	if(mxGetField( prhs[0], 0, "yr" )!=NULL)
		{
		set_yr = true;
		yr = mxGetPr( mxGetField( prhs[0], 0, "yr" ) );
		}
	// yr_e
	if(mxGetField( prhs[0], 0, "yr_e" )!=NULL)
		{
		set_yr_e = true;
		yr_e = mxGetPr( mxGetField( prhs[0], 0, "yr_e" ) );
		}
	// Z
	if(mxGetField( prhs[0], 0, "Z" )!=NULL)
		{
		set_Z = true;
		Z = mxGetPr( mxGetField( prhs[0], 0, "Z" ) );
		}
	// Z_e
	if(mxGetField( prhs[0], 0, "Z_e" )!=NULL)
		{
		set_Z_e = true;
		Z_e = mxGetPr( mxGetField( prhs[0], 0, "Z_e" ) );
		}
	// Zl
	if(mxGetField( prhs[0], 0, "Zl" )!=NULL)
		{
		set_Zl = true;
		Zl = mxGetPr( mxGetField( prhs[0], 0, "Zl" ) );
		}
	// Zl_e
	if(mxGetField( prhs[0], 0, "Zl_e" )!=NULL)
		{
		set_Zl_e = true;
		Zl_e = mxGetPr( mxGetField( prhs[0], 0, "Zl_e" ) );
		}
	// Zu
	if(mxGetField( prhs[0], 0, "Zu" )!=NULL)
		{
		set_Zu = true;
		Zu = mxGetPr( mxGetField( prhs[0], 0, "Zu" ) );
		}
	// Zu_e
	if(mxGetField( prhs[0], 0, "Zu_e" )!=NULL)
		{
		set_Zu_e = true;
		Zu_e = mxGetPr( mxGetField( prhs[0], 0, "Zu_e" ) );
		}
	// z
	if(mxGetField( prhs[0], 0, "z" )!=NULL)
		{
		set_z = true;
		z = mxGetPr( mxGetField( prhs[0], 0, "z" ) );
		}
	// z_e
	if(mxGetField( prhs[0], 0, "z_e" )!=NULL)
		{
		set_z_e = true;
		z_e = mxGetPr( mxGetField( prhs[0], 0, "z_e" ) );
		}
	// zl
	if(mxGetField( prhs[0], 0, "zl" )!=NULL)
		{
		set_zl = true;
		zl = mxGetPr( mxGetField( prhs[0], 0, "zl" ) );
		}
	// zl_e
	if(mxGetField( prhs[0], 0, "zl_e" )!=NULL)
		{
		set_zl_e = true;
		zl_e = mxGetPr( mxGetField( prhs[0], 0, "zl_e" ) );
		}
	// zu
	if(mxGetField( prhs[0], 0, "zu" )!=NULL)
		{
		set_zu = true;
		zu = mxGetPr( mxGetField( prhs[0], 0, "zu" ) );
		}
	// zu_e
	if(mxGetField( prhs[0], 0, "zu_e" )!=NULL)
		{
		set_zu_e = true;
		zu_e = mxGetPr( mxGetField( prhs[0], 0, "zu_e" ) );
		}
	// constr
	// constr type
	if(mxGetField( prhs[0], 0, "constr_type" )!=NULL)
		{
		constr_type = mxArrayToString( mxGetField( prhs[0], 0, "constr_type" ) );
		}
	// x0
	if(mxGetField( prhs[0], 0, "x0" )!=NULL)
		{
		set_x0 = true;
		x0 = mxGetPr( mxGetField( prhs[0], 0, "x0" ) );
		}
	// Jbx
	if(mxGetField( prhs[0], 0, "Jbx" )!=NULL)
		{
		set_Jbx = true;
		Jbx = mxGetPr( mxGetField( prhs[0], 0, "Jbx" ) );
		}
	// lbx
	if(mxGetField( prhs[0], 0, "lbx" )!=NULL)
		{
		set_lbx = true;
		lbx = mxGetPr( mxGetField( prhs[0], 0, "lbx" ) );
		}
	// ubx
	if(mxGetField( prhs[0], 0, "ubx" )!=NULL)
		{
		set_ubx = true;
		ubx = mxGetPr( mxGetField( prhs[0], 0, "ubx" ) );
		}
	// Jbu
	if(mxGetField( prhs[0], 0, "Jbu" )!=NULL)
		{
		set_Jbu = true;
		Jbu = mxGetPr( mxGetField( prhs[0], 0, "Jbu" ) );
		}
	// lbu
	if(mxGetField( prhs[0], 0, "lbu" )!=NULL)
		{
		set_lbu = true;
		lbu = mxGetPr( mxGetField( prhs[0], 0, "lbu" ) );
		}
	// ubu
	if(mxGetField( prhs[0], 0, "ubu" )!=NULL)
		{
		set_ubu = true;
		ubu = mxGetPr( mxGetField( prhs[0], 0, "ubu" ) );
		}
	// C
	if(mxGetField( prhs[0], 0, "C" )!=NULL)
		{
		set_C = true;
		C = mxGetPr( mxGetField( prhs[0], 0, "C" ) );
		}
	// D
	if(mxGetField( prhs[0], 0, "D" )!=NULL)
		{
		set_D = true;
		D = mxGetPr( mxGetField( prhs[0], 0, "D" ) );
		}
	// lg
	if(mxGetField( prhs[0], 0, "lg" )!=NULL)
		{
		set_lg = true;
		lg = mxGetPr( mxGetField( prhs[0], 0, "lg" ) );
		}
	// ug
	if(mxGetField( prhs[0], 0, "ug" )!=NULL)
		{
		set_ug = true;
		ug = mxGetPr( mxGetField( prhs[0], 0, "ug" ) );
		}
	// C_e
	if(mxGetField( prhs[0], 0, "C_e" )!=NULL)
		{
		set_C_e = true;
		C_e = mxGetPr( mxGetField( prhs[0], 0, "C_e" ) );
		}
	// lg_e
	if(mxGetField( prhs[0], 0, "lg_e" )!=NULL)
		{
		set_lg_e = true;
		lg_e = mxGetPr( mxGetField( prhs[0], 0, "lg_e" ) );
		}
	// ug_e
	if(mxGetField( prhs[0], 0, "ug_e" )!=NULL)
		{
		set_ug_e = true;
		ug_e = mxGetPr( mxGetField( prhs[0], 0, "ug_e" ) );
		}
	// lh
	if(mxGetField( prhs[0], 0, "lh" )!=NULL)
		{
		set_lh = true;
		lh = mxGetPr( mxGetField( prhs[0], 0, "lh" ) );
		}
	// uh
	if(mxGetField( prhs[0], 0, "uh" )!=NULL)
		{
		set_uh = true;
		uh = mxGetPr( mxGetField( prhs[0], 0, "uh" ) );
		}
	// lh_e
	if(mxGetField( prhs[0], 0, "lh_e" )!=NULL)
		{
		set_lh_e = true;
		lh_e = mxGetPr( mxGetField( prhs[0], 0, "lh_e" ) );
		}
	// uh_e
	if(mxGetField( prhs[0], 0, "uh_e" )!=NULL)
		{
		set_uh_e = true;
		uh_e = mxGetPr( mxGetField( prhs[0], 0, "uh_e" ) );
		}
	// Jsbu
	if(mxGetField( prhs[0], 0, "Jsbu" )!=NULL)
		{
		set_Jsbu = true;
		Jsbu = mxGetPr( mxGetField( prhs[0], 0, "Jsbu" ) );
		}
	// lsbu
//	if(mxGetField( prhs[0], 0, "lsbu" )!=NULL)
//		{
//		set_lsbu = true;
//		lsbu = mxGetPr( mxGetField( prhs[0], 0, "lsbu" ) );
//		}
	// usbu
//	if(mxGetField( prhs[0], 0, "usbu" )!=NULL)
//		{
//		set_usbu = true;
//		usbu = mxGetPr( mxGetField( prhs[0], 0, "usbu" ) );
//		}
	// Jsbx
	if(mxGetField( prhs[0], 0, "Jsbx" )!=NULL)
		{
		set_Jsbx = true;
		Jsbx = mxGetPr( mxGetField( prhs[0], 0, "Jsbx" ) );
		}
	// lsbx
//	if(mxGetField( prhs[0], 0, "lsbx" )!=NULL)
//		{
//		set_lsbx = true;
//		lsbx = mxGetPr( mxGetField( prhs[0], 0, "lsbx" ) );
//		}
	// usbx
//	if(mxGetField( prhs[0], 0, "usbx" )!=NULL)
//		{
//		set_usbx = true;
//		usbx = mxGetPr( mxGetField( prhs[0], 0, "usbx" ) );
//		}
	// Jsg
	if(mxGetField( prhs[0], 0, "Jsg" )!=NULL)
		{
		set_Jsg = true;
		Jsg = mxGetPr( mxGetField( prhs[0], 0, "Jsg" ) );
		}
	// lsg
//	if(mxGetField( prhs[0], 0, "lsg" )!=NULL)
//		{
//		set_lsg = true;
//		lsg = mxGetPr( mxGetField( prhs[0], 0, "lsg" ) );
//		}
	// usg
//	if(mxGetField( prhs[0], 0, "usg" )!=NULL)
//		{
//		set_usg = true;
//		usg = mxGetPr( mxGetField( prhs[0], 0, "usg" ) );
//		}
	// Jsg_e
	if(mxGetField( prhs[0], 0, "Jsg_e" )!=NULL)
		{
		set_Jsg_e = true;
		Jsg_e = mxGetPr( mxGetField( prhs[0], 0, "Jsg_e" ) );
		}
	// lsg_e
//	if(mxGetField( prhs[0], 0, "lsg_e" )!=NULL)
//		{
//		set_lsg_e = true;
//		lsg_e = mxGetPr( mxGetField( prhs[0], 0, "lsg_e" ) );
//		}
	// usg_e
//	if(mxGetField( prhs[0], 0, "usg_e" )!=NULL)
//		{
//		set_usg_e = true;
//		usg_e = mxGetPr( mxGetField( prhs[0], 0, "usg_e" ) );
//		}
	// Jsh
	if(mxGetField( prhs[0], 0, "Jsh" )!=NULL)
		{
		set_Jsh = true;
		Jsh = mxGetPr( mxGetField( prhs[0], 0, "Jsh" ) );
		}
	// lsh
//	if(mxGetField( prhs[0], 0, "lsh" )!=NULL)
//		{
//		set_lsh = true;
//		lsh = mxGetPr( mxGetField( prhs[0], 0, "lsh" ) );
//		}
	// ush
//	if(mxGetField( prhs[0], 0, "ush" )!=NULL)
//		{
//		set_ush = true;
//		ush = mxGetPr( mxGetField( prhs[0], 0, "ush" ) );
//		}
	// Jsh_e
	if(mxGetField( prhs[0], 0, "Jsh_e" )!=NULL)
		{
		set_Jsh_e = true;
		Jsh_e = mxGetPr( mxGetField( prhs[0], 0, "Jsh_e" ) );
		}
	// lsh_e
//	if(mxGetField( prhs[0], 0, "lsh_e" )!=NULL)
//		{
//		set_lsh_e = true;
//		lsh_e = mxGetPr( mxGetField( prhs[0], 0, "lsh_e" ) );
//		}
	// ush_e
//	if(mxGetField( prhs[0], 0, "ush_e" )!=NULL)
//		{
//		set_ush_e = true;
//		ush_e = mxGetPr( mxGetField( prhs[0], 0, "ush_e" ) );
//		}
	// dyn
	// dyn type
	if(mxGetField( prhs[0], 0, "dyn_type" )!=NULL)
		{
		dyn_type = mxArrayToString( mxGetField( prhs[0], 0, "dyn_type" ) );
		}
	// trajectory initialization
	// u_init
	if(mxGetField( prhs[0], 0, "u_init" )!=NULL)
		{
		set_u_init = true;
		u_init = mxGetPr( mxGetField( prhs[0], 0, "u_init" ) );
		}
	// x_init
	if(mxGetField( prhs[0], 0, "x_init" )!=NULL)
		{
		set_x_init = true;
		x_init = mxGetPr( mxGetField( prhs[0], 0, "x_init" ) );
		}


	// opts_struct

	int N;							bool set_param_scheme_N = false;
	char * nlp_solver;
	char * qp_solver;
	int qp_solver_N_pcond;			bool set_qp_solver_N_pcond = false;
	char *sim_method;
	int sim_method_num_stages;		bool set_sim_method_num_stages = false;
	int sim_method_num_steps;		bool set_sim_method_num_steps = false;

	// param_scheme_NN
	if(mxGetField( prhs[1], 0, "param_scheme_N" )!=NULL)
		{
		set_param_scheme_N = true;
		N = mxGetScalar( mxGetField( prhs[1], 0, "param_scheme_N" ) );
		}
	else
		{
		mexPrintf("\nerror: ocp_create: param_scheme_N not set!\n");
		return;
		}
	// nlp_solver
	// TODO check
	nlp_solver = mxArrayToString( mxGetField( prhs[1], 0, "nlp_solver" ) );
	// qp_solver
	// TODO check
	qp_solver = mxArrayToString( mxGetField( prhs[1], 0, "qp_solver" ) );
	// N_part_cond
	if(mxGetField( prhs[1], 0, "qp_solver_N_pcond" )!=NULL)
		{
		set_qp_solver_N_pcond = true;
		qp_solver_N_pcond = mxGetScalar( mxGetField( prhs[1], 0, "qp_solver_N_pcond" ) );
		}
	// sim_method
	// TODO check
	sim_method = mxArrayToString( mxGetField( prhs[1], 0, "sim_method" ) );
	// sim_method_num_stages
	if(mxGetField( prhs[1], 0, "sim_method_num_stages" )!=NULL)
		{
		set_sim_method_num_stages = true;
		sim_method_num_stages = mxGetScalar( mxGetField( prhs[1], 0, "sim_method_num_stages" ) );
		}
	// sim_method_num_steps
	if(mxGetField( prhs[1], 0, "sim_method_num_steps" )!=NULL)
		{
		set_sim_method_num_steps = true;
		sim_method_num_steps = mxGetScalar( mxGetField( prhs[1], 0, "sim_method_num_steps" ) );
		}



	/* LHS */

	// field names of output struct
	char *fieldnames[6];
	fieldnames[0] = (char*) mxMalloc(50);
	fieldnames[1] = (char*) mxMalloc(50);
	fieldnames[2] = (char*) mxMalloc(50);
	fieldnames[3] = (char*) mxMalloc(50);
	fieldnames[4] = (char*) mxMalloc(50);
	fieldnames[5] = (char*) mxMalloc(50);

	memcpy(fieldnames[0],"config",sizeof("config"));
	memcpy(fieldnames[1],"dims",sizeof("dims"));
	memcpy(fieldnames[2],"opts",sizeof("opts"));
	memcpy(fieldnames[3],"in",sizeof("in"));
	memcpy(fieldnames[4],"out",sizeof("out"));
	memcpy(fieldnames[5],"solver",sizeof("solver"));

	// create output struct
	plhs[0] = mxCreateStructMatrix(1, 1, 6, (const char **) fieldnames);

	mxFree( fieldnames[0] );
	mxFree( fieldnames[1] );
	mxFree( fieldnames[2] );
	mxFree( fieldnames[3] );
	mxFree( fieldnames[4] );
	mxFree( fieldnames[5] );



	/* plan & config */

	ocp_nlp_plan *plan = ocp_nlp_plan_create(N);

	// nlp solver
	if(!strcmp(nlp_solver, "sqp"))
		{
		plan->nlp_solver = SQP;
		}
	else if(!strcmp(nlp_solver, "sqp_rti"))
		{
		plan->nlp_solver = SQP_RTI;
		}
	else
		{
		mexPrintf("\nnlp_solver not supported: %s\n", nlp_solver);
		return;
		}
	
	// cost
	if(!strcmp(cost_type, "linear_ls"))
		{
		for(ii=0; ii<N; ii++)
			{
			plan->nlp_cost[ii] = LINEAR_LS;
			}
		}
	else if(!strcmp(cost_type, "nonlinear_ls"))
		{
		for(ii=0; ii<N; ii++)
			{
			plan->nlp_cost[ii] = NONLINEAR_LS;
			}
		}
	else if(!strcmp(cost_type, "ext_cost"))
		{
		for(ii=0; ii<N; ii++)
			{
			plan->nlp_cost[ii] = EXTERNALLY_PROVIDED;
			}
		}
	else
		{
		mexPrintf("\ncost_type not supported: %s\n", cost_type);
		return;
		}
	if(!strcmp(cost_e_type, "linear_ls"))
		{
		plan->nlp_cost[N] = LINEAR_LS;
		}
	else if(!strcmp(cost_e_type, "nonlinear_ls"))
		{
		plan->nlp_cost[N] = NONLINEAR_LS;
		}
	else if(!strcmp(cost_e_type, "ext_cost"))
		{
		plan->nlp_cost[N] = EXTERNALLY_PROVIDED;
		}
	else
		{
		mexPrintf("\ncost_e_type not supported: %s\n", cost_e_type);
		return;
		}

	// dynamics
	if(!strcmp(dyn_type, "explicit"))
		{
		for(ii=0; ii<N; ii++)
			{
			plan->nlp_dynamics[ii] = CONTINUOUS_MODEL;
			}
		if(!strcmp(sim_method, "erk"))
			{
			for(ii=0; ii<N; ii++)
				{
				plan->sim_solver_plan[ii].sim_solver = ERK;
				}
			}
		else
			{
			mexPrintf("\nsim_method not supported for explicit dynamics: %s\n", sim_method);
			return;
			}
		}
	else if(!strcmp(dyn_type, "implicit"))
		{
		for(ii=0; ii<N; ii++)
			{
			plan->nlp_dynamics[ii] = CONTINUOUS_MODEL;
			}
		if(!strcmp(sim_method, "irk"))
			{
			for(ii=0; ii<N; ii++)
				{
				plan->sim_solver_plan[ii].sim_solver = IRK;
				}
			}
		else
			{
			mexPrintf("\nsim_method not supported for implicit dynamics: %s\n", sim_method);
			return;
			}
//		plan->sim_solver_plan[ii].sim_solver = LIFTED_IRK;
		}
	else // TODO gnsf / discrete / linear
		{
//		plan->nlp_dynamics[ii] = DISCRETE_MODEL;
//		plan->sim_solver_plan[ii].sim_solver = GNSF;
		mexPrintf("\ndyn_type not supported %s\n", dyn_type);
		return;
		}
	
	// constraints
	if(!strcmp(constr_type, "bgh"))
		{
		for(ii=0; ii<=N; ii++)
			{
			plan->nlp_constraints[ii] = BGH;
			}
		}
	else
		{
		mexPrintf("\nconstr_type not supported: %s\n", constr_type);
		return;
		}

	// qp solver
	if(!strcmp(qp_solver, "partial_condensing_hpipm"))
		{
		plan->ocp_qp_solver_plan.qp_solver = PARTIAL_CONDENSING_HPIPM;
		}
	else if(!strcmp(qp_solver, "full_condensing_hpipm"))
		{
		plan->ocp_qp_solver_plan.qp_solver = FULL_CONDENSING_HPIPM;
		}
	else
		{
		mexPrintf("\nqp_solver not supported: %s\n", qp_solver);
		return;
		}

	// TODO checks on initialization of plan !!!!!!!!

    ocp_nlp_config *config = ocp_nlp_config_create(*plan);


	ocp_nlp_plan_destroy(plan);
	


	/* dims */

	ocp_nlp_dims *dims = ocp_nlp_dims_create(config);
	// allocate tmp
	i_ptr = (int *) malloc((N+1)*sizeof(int));
	// nx
	for(ii=0; ii<=N; ii++)
		i_ptr[ii] = nx;
	ocp_nlp_dims_set_opt_vars(config, dims, "nx", i_ptr);
	// nu
	for(ii=0; ii<N; ii++)
		i_ptr[ii] = nu;
	i_ptr[N] = 0;
	ocp_nlp_dims_set_opt_vars(config, dims, "nu", i_ptr);
	// ns
	if(ns!=nsbu+nsbx+nsg+nsh)
		{
		mexPrintf("\nerror: ns!=nsbu+nsbx+nsg+nsh\n");
		return;
		}
	if(ns_e!=nsbx+nsg_e+nsh_e)
		{
		mexPrintf("\nerror: ns_e!=nsbx+nsg_e+nsh_e\n");
		return;
		}
	// TODO fix stage 0 !!!!!!!!
	i_ptr[0] = nsbu+nsg+nsh; // XXX not nsbx !!!!!!!!!!
	for(ii=1; ii<N; ii++)
		i_ptr[ii] = ns;
	i_ptr[N] = ns_e;
	ocp_nlp_dims_set_opt_vars(config, dims, "ns", i_ptr);
	// free tmp
	free(i_ptr);
	// ny
	if(set_ny)
		{
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_dims_set_cost(config, dims, ii, "ny", &ny);
			}
		}
	// ny_e
	if(set_ny_e)
		{
		ocp_nlp_dims_set_cost(config, dims, N, "ny", &ny_e);
		}
	// nbx
	ocp_nlp_dims_set_constraints(config, dims, 0, "nbx", &nx);
	if(set_nbx)
		{
		for(ii=1; ii<=N; ii++)
			{
			ocp_nlp_dims_set_constraints(config, dims, ii, "nbx", &nbx);
			}
		}
	// nbu
	if(set_nbu)
		{
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_dims_set_constraints(config, dims, ii, "nbu", &nbu);
			}
		}
	// ng
	if(set_ng)
		{
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_dims_set_constraints(config, dims, ii, "ng", &ng);
			}
		}
	// ng_e
	if(set_ng_e)
		{
		ocp_nlp_dims_set_constraints(config, dims, N, "ng", &ng_e);
		}
	// nh
	if(set_nh)
		{
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_dims_set_constraints(config, dims, ii, "nh", &nh);
			}
		}
	// nh_e
	if(set_nh_e)
		{
		ocp_nlp_dims_set_constraints(config, dims, N, "nh", &nh_e);
		}
	// nsbx
	if(set_nsbx)
		{
//		for(ii=0; ii<=N; ii++)
		for(ii=1; ii<=N; ii++) // TODO fix stage 0 !!!!!
			{
			ocp_nlp_dims_set_constraints(config, dims, ii, "nsbx", &nsbx);
			}
		}
	// nsbu
	if(set_nsbu)
		{
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_dims_set_constraints(config, dims, ii, "nsbu", &nsbu);
			}
		}
	// nsg
	if(set_nsg)
		{
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_dims_set_constraints(config, dims, ii, "nsg", &nsg);
			}
		}
	// nsg_e
	if(set_nsg_e)
		{
		ocp_nlp_dims_set_constraints(config, dims, N, "nsg", &nsg_e);
		}
	// nsh
	if(set_nsh)
		{
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_dims_set_constraints(config, dims, ii, "nsh", &nsh);
			}
		}
	// nsh_e
	if(set_nsh_e)
		{
		ocp_nlp_dims_set_constraints(config, dims, N, "nsh", &nsh_e);
		}
			


	/* opts */

	void *opts = ocp_nlp_opts_create(config, dims);

	// qp_solver_N_pcond
	if(set_qp_solver_N_pcond)
		{
		ocp_nlp_opts_set(config, opts, "pcond_N2", &qp_solver_N_pcond);
		}
	// sim_method_num_stages
	if(set_sim_method_num_stages)
		{
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_dynamics_opts_set(config, opts, ii, "num_stages", &sim_method_num_stages);
			}
		}
	// sim_method_num_steps
	if(set_sim_method_num_steps)
		{
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_dynamics_opts_set(config, opts, ii, "num_steps", &sim_method_num_steps);
			}
		}




	/* in */

	ocp_nlp_in *in = ocp_nlp_in_create(config, dims);

	// shooting nodes
	double Ts = T/N;
	for(ii=0; ii<N; ii++)
		{
		ocp_nlp_in_set(config, dims, in, ii, "Ts", &Ts);
		ocp_nlp_cost_model_set(config, dims, in, ii, "scaling", &Ts);
		}

	// cost: ls
	if(!strcmp(cost_type, "linear_ls"))
		{
		// lagrange term
		if(set_Vu)
			{
			for(ii=0; ii<N; ii++)
				{
				ocp_nlp_cost_model_set(config, dims, in, ii, "Vu", Vu);
				}
			}
		if(set_Vx)
			{
			for(ii=0; ii<N; ii++)
				{
				ocp_nlp_cost_model_set(config, dims, in, ii, "Vx", Vx);
				}
			}
		if(set_W)
			{
			for(ii=0; ii<N; ii++)
				{
				ocp_nlp_cost_model_set(config, dims, in, ii, "W", W);
				}
			}
		if(set_yr)
			{
			for(ii=0; ii<N; ii++)
				{
				ocp_nlp_cost_model_set(config, dims, in, ii, "y_ref", yr);
				}
			}
		}
	if(!strcmp(cost_e_type, "linear_ls"))
		{
		// mayer term
		if(set_Vx_e)
			{
			ocp_nlp_cost_model_set(config, dims, in, N, "Vx", Vx_e);
			}
		if(set_W_e)
			{
			ocp_nlp_cost_model_set(config, dims, in, N, "W", W_e);
			}
		if(set_yr_e)
			{
			ocp_nlp_cost_model_set(config, dims, in, N, "y_ref", yr_e);
			}
		}
	// cost: nls
	if(!strcmp(cost_type, "nonlinear_ls"))
		{
		// lagrange term
		if(set_W)
			{
			for(ii=0; ii<N; ii++)
				{
				ocp_nlp_cost_model_set(config, dims, in, ii, "W", W);
				}
			}
		if(set_yr)
			{
			for(ii=0; ii<N; ii++)
				{
				ocp_nlp_cost_model_set(config, dims, in, ii, "y_ref", yr);
				}
			}
		}
	if(!strcmp(cost_e_type, "nonlinear_ls"))
		{
		// mayer term
		if(set_W_e)
			{
			ocp_nlp_cost_model_set(config, dims, in, N, "W", W_e);
			}
		if(set_yr_e)
			{
			ocp_nlp_cost_model_set(config, dims, in, N, "y_ref", yr_e);
			}
		}
	// slacks
	if(set_Z)
		{
		d_ptr = malloc(ns*sizeof(double));
		for(ii=0; ii<ns; ii++)
			{
			d_ptr[ii] = Z[ii+ns*ii];
			}
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_cost_model_set(config, dims, in, ii, "Z", d_ptr);
			}
		free(d_ptr);
		}
	if(set_Z_e)
		{
		d_ptr = malloc(ns*sizeof(double));
		for(ii=0; ii<ns; ii++)
			{
			d_ptr[ii] = Z_e[ii+ns*ii];
			}
		ocp_nlp_cost_model_set(config, dims, in, N, "Z", d_ptr);
		free(d_ptr);
		}
	if(set_Zl)
		{
		d_ptr = malloc(ns*sizeof(double));
		for(ii=0; ii<ns; ii++)
			{
			d_ptr[ii] = Zl[ii+ns*ii];
			}
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_cost_model_set(config, dims, in, ii, "Zl", d_ptr);
			}
		free(d_ptr);
		}
	if(set_Zl_e)
		{
		d_ptr = malloc(ns*sizeof(double));
		for(ii=0; ii<ns; ii++)
			{
			d_ptr[ii] = Zl_e[ii+ns*ii];
			}
		ocp_nlp_cost_model_set(config, dims, in, N, "Zl", d_ptr);
		free(d_ptr);
		}
	if(set_Zu)
		{
		d_ptr = malloc(ns*sizeof(double));
		for(ii=0; ii<ns; ii++)
			{
			d_ptr[ii] = Zu[ii+ns*ii];
			}
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_cost_model_set(config, dims, in, ii, "Zu", d_ptr);
			}
		free(d_ptr);
		}
	if(set_Zu_e)
		{
		d_ptr = malloc(ns*sizeof(double));
		for(ii=0; ii<ns; ii++)
			{
			d_ptr[ii] = Zu_e[ii+ns*ii];
			}
		ocp_nlp_cost_model_set(config, dims, in, N, "Zu", d_ptr);
		free(d_ptr);
		}
	if(set_z)
		{
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_cost_model_set(config, dims, in, ii, "z", z);
			}
		}
	if(set_z_e)
		{
		ocp_nlp_cost_model_set(config, dims, in, N, "z", z_e);
		}
	if(set_zl)
		{
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_cost_model_set(config, dims, in, ii, "zl", zl);
			}
		}
	if(set_zl_e)
		{
		ocp_nlp_cost_model_set(config, dims, in, N, "zl", zl_e);
		}
	if(set_zu)
		{
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_cost_model_set(config, dims, in, ii, "zu", zu);
			}
		}
	if(set_zu_e)
		{
		ocp_nlp_cost_model_set(config, dims, in, N, "zu", zu_e);
		}

	// constraints: bgh

	double acados_inf = 1e8;

	// x0 is always bounded on all components !!!
	i_ptr = malloc(nx*sizeof(int));
	for(ii=0; ii<nx; ii++)
		{
		i_ptr[ii] = ii;
		}
	ocp_nlp_constraints_model_set(config, dims, in, 0, "idxbx", i_ptr);
	free(i_ptr);
	if(set_x0)
		{
		ocp_nlp_constraints_model_set(config, dims, in, 0, "lbx", x0);
		ocp_nlp_constraints_model_set(config, dims, in, 0, "ubx", x0);
		}
	else
		{
		d_ptr = malloc(nx*sizeof(double));
		for(ii=0; ii<nx; ii++)
			{
			d_ptr[ii] = - acados_inf;
			}
		ocp_nlp_constraints_model_set(config, dims, in, 0, "lbx", d_ptr);
		for(ii=0; ii<nx; ii++)
			{
			d_ptr[ii] = acados_inf;
			}
		ocp_nlp_constraints_model_set(config, dims, in, 0, "ubx", d_ptr);
		free(d_ptr);
		}
	tmp_idx = malloc(nbx*sizeof(int));
	if(set_Jbx)
		{
		for(ii=0; ii<nbx; ii++)
			{
			idx = -1;
			for(jj=0; jj<nx; jj++)
				{
				if(Jbx[ii+nbx*jj]!=0.0)
					{
					tmp_idx[ii] = jj;
					idx = jj;
					}
				}
			}
		ii = 1;
		for(; ii<=N; ii++)
			{
			ocp_nlp_constraints_model_set(config, dims, in, ii, "idxbx", tmp_idx);
			}
		}
	if(set_lbx)
		{
		if(!set_x0)
			{
			d_ptr = malloc(nx*sizeof(double));
			for(ii=0; ii<nx; ii++)
				{
				d_ptr[ii] = - acados_inf;
				}
			for(ii=0; ii<nbx; ii++)
				{
				d_ptr[tmp_idx[ii]] = lbx[ii];
				}
			ocp_nlp_constraints_model_set(config, dims, in, 0, "lbx", d_ptr);
			free(d_ptr);
			}
		ii = 1;
		for(; ii<=N; ii++)
			{
			ocp_nlp_constraints_model_set(config, dims, in, ii, "lbx", lbx);
			}
		}
	if(set_ubx)
		{
		if(!set_x0)
			{
			d_ptr = malloc(nx*sizeof(double));
			for(ii=0; ii<nx; ii++)
				{
				d_ptr[ii] = acados_inf;
				}
			for(ii=0; ii<nbx; ii++)
				{
				d_ptr[tmp_idx[ii]] = ubx[ii];
				}
			ocp_nlp_constraints_model_set(config, dims, in, 0, "ubx", d_ptr);
			free(d_ptr);
			}
		ii = 1;
		for(; ii<=N; ii++)
			{
			ocp_nlp_constraints_model_set(config, dims, in, ii, "ubx", ubx);
			}
		}
	free(tmp_idx);
	if(set_Jbu)
		{
		i_ptr = malloc(nbu*sizeof(int));
		for(ii=0; ii<nbu; ii++)
			{
			idx = -1;
			for(jj=0; jj<nu; jj++)
				{
				if(Jbu[ii+nbu*jj]!=0.0)
					{
					i_ptr[ii] = jj;
					idx = jj;
					}
				}
			}
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_constraints_model_set(config, dims, in, ii, "idxbu", i_ptr);
			}
		free(i_ptr);
		}
	if(set_lbu)
		{
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_constraints_model_set(config, dims, in, ii, "lbu", lbu);
			}
		}
	if(set_ubu)
		{
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_constraints_model_set(config, dims, in, ii, "ubu", ubu);
			}
		}
	if(set_C)
		{
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_constraints_model_set(config, dims, in, ii, "C", C);
			}
		}
	if(set_D)
		{
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_constraints_model_set(config, dims, in, ii, "D", D);
			}
		}
	if(set_lg)
		{
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_constraints_model_set(config, dims, in, ii, "lg", lg);
			}
		}
	if(set_ug)
		{
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_constraints_model_set(config, dims, in, ii, "ug", ug);
			}
		}
	if(set_C_e)
		{
		ocp_nlp_constraints_model_set(config, dims, in, N, "C", C_e);
		}
	if(set_lg_e)
		{
		ocp_nlp_constraints_model_set(config, dims, in, N, "lg", lg_e);
		}
	if(set_ug_e)
		{
		ocp_nlp_constraints_model_set(config, dims, in, N, "ug", ug_e);
		}
	if(set_lh)
		{
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_constraints_model_set(config, dims, in, ii, "lh", lh);
			}
		}
	if(set_lh_e)
		{
		ocp_nlp_constraints_model_set(config, dims, in, N, "lh", lh_e);
		}
	if(set_uh)
		{
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_constraints_model_set(config, dims, in, ii, "uh", uh);
			}
		}
	if(set_uh_e)
		{
		ocp_nlp_constraints_model_set(config, dims, in, N, "uh", uh_e);
		}
	if(set_Jsbu)
		{
		i_ptr = malloc(nsbu*sizeof(int));
		for(ii=0; ii<nsbu; ii++)
			{
			idx = -1;
			for(jj=0; jj<nbu; jj++)
				{
				if(Jsbu[jj+nbu*ii]!=0.0)
					{
					i_ptr[ii] = jj;
					idx = jj;
					}
				}
			}
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_constraints_model_set(config, dims, in, ii, "idxsbu", i_ptr);
			}
		free(i_ptr);
		}
//	if(set_lsbu)
//		{
//		for(ii=0; ii<N; ii++)
//			{
//			ocp_nlp_constraints_model_set(config, dims, in, ii, "lsbu", lsbu);
//			}
//		}
//	if(set_usbu)
//		{
//		for(ii=0; ii<N; ii++)
//			{
//			ocp_nlp_constraints_model_set(config, dims, in, ii, "usbu", usbu);
//			}
//		}
	if(set_Jsbx)
		{
		i_ptr = malloc(nsbx*sizeof(int));
		for(ii=0; ii<nsbx; ii++)
			{
			idx = -1;
			for(jj=0; jj<nbx; jj++)
				{
				if(Jsbx[jj+nbx*ii]!=0.0)
					{
					i_ptr[ii] = jj;
					idx = jj;
					}
				}
			}
//		for(ii=0; ii<=N; ii++)
		for(ii=1; ii<=N; ii++) // TODO stage 0 !!!!!!!!!!
			{
			ocp_nlp_constraints_model_set(config, dims, in, ii, "idxsbx", i_ptr);
			}
		free(i_ptr);
		}
//	if(set_lsbx)
//		{
//		for(ii=0; ii<=N; ii++)
//		for(ii=1; ii<=N; ii++) // TODO stage 0 !!!!!!!!!!
//			{
//			ocp_nlp_constraints_model_set(config, dims, in, ii, "lsbx", lsbx);
//			}
//		}
//	if(set_usbx)
//		{
//		for(ii=0; ii<=N; ii++)
//		for(ii=1; ii<=N; ii++) // TODO stage 0 !!!!!!!!!!
//			{
//			ocp_nlp_constraints_model_set(config, dims, in, ii, "usbx", usbx);
//			}
//		}
	if(set_Jsg)
		{
		i_ptr = malloc(nsg*sizeof(int));
		for(ii=0; ii<nsg; ii++)
			{
			idx = -1;
			for(jj=0; jj<ng; jj++)
				{
				if(Jsg[jj+ng*ii]!=0.0)
					{
					i_ptr[ii] = jj;
					idx = jj;
					}
				}
			}
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_constraints_model_set(config, dims, in, ii, "idxsg", i_ptr);
			}
		free(i_ptr);
		}
//	if(set_lsg)
//		{
//		for(ii=0; ii<N; ii++)
//			{
//			ocp_nlp_constraints_model_set(config, dims, in, ii, "lsg", lsg);
//			}
//		}
//	if(set_usg)
//		{
//		for(ii=0; ii<N; ii++)
//			{
//			ocp_nlp_constraints_model_set(config, dims, in, ii, "usg", usg);
//			}
//		}
	if(set_Jsg_e)
		{
		i_ptr = malloc(nsg_e*sizeof(int));
		for(ii=0; ii<nsg_e; ii++)
			{
			idx = -1;
			for(jj=0; jj<ng_e; jj++)
				{
				if(Jsg_e[jj+ng_e*ii]!=0.0)
					{
					i_ptr[ii] = jj;
					idx = jj;
					}
				}
			}
		ocp_nlp_constraints_model_set(config, dims, in, N, "idxsg", i_ptr);
		free(i_ptr);
		}
//	if(set_lsg_e)
//		{
//		ocp_nlp_constraints_model_set(config, dims, in, N, "lsg", lsg_e);
//		}
//	if(set_usg_e)
//		{
//		ocp_nlp_constraints_model_set(config, dims, in, N, "usg", usg_e);
//		}
	if(set_Jsh)
		{
		i_ptr = malloc(nsh*sizeof(int));
		for(ii=0; ii<nsh; ii++)
			{
			idx = -1;
			for(jj=0; jj<nh; jj++)
				{
				if(Jsh[jj+nh*ii]!=0.0)
					{
					i_ptr[ii] = jj;
					idx = jj;
					}
				}
			}
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_constraints_model_set(config, dims, in, ii, "idxsh", i_ptr);
			}
		free(i_ptr);
		}
//	if(set_lsh)
//		{
//		for(ii=0; ii<N; ii++)
//			{
//			ocp_nlp_constraints_model_set(config, dims, in, ii, "lsh", lsh);
//			}
//		}
//	if(set_ush)
//		{
//		for(ii=0; ii<N; ii++)
//			{
//			ocp_nlp_constraints_model_set(config, dims, in, ii, "ush", ush);
//			}
//		}
	if(set_Jsh_e)
		{
		i_ptr = malloc(nsh_e*sizeof(int));
		for(ii=0; ii<nsh_e; ii++)
			{
			idx = -1;
			for(jj=0; jj<nh_e; jj++)
				{
				if(Jsh_e[jj+nh_e*ii]!=0.0)
					{
					i_ptr[ii] = jj;
					idx = jj;
					}
				}
			}
		ocp_nlp_constraints_model_set(config, dims, in, N, "idxsh", i_ptr);
		free(i_ptr);
		}
//	if(set_lsh_e)
//		{
//		ocp_nlp_constraints_model_set(config, dims, in, N, "lsh", lsh_e);
//		}
//	if(set_ush_e)
//		{
//		ocp_nlp_constraints_model_set(config, dims, in, N, "ush", ush_e);
//		}



	/* out */

	ocp_nlp_out *out = ocp_nlp_out_create(config, dims);

	if(set_x_init)
		{
		for(ii=0; ii<=N; ii++)
			{
			ocp_nlp_out_set(config, dims, out, ii, "x", x_init+ii*nx);
			}
		}
	else // initialize to zero
		{
		double *x_init = calloc(nx, sizeof(double));
		for(ii=0; ii<=N; ii++)
			{
			ocp_nlp_out_set(config, dims, out, ii, "x", x_init);
			}
		free(x_init);
		}
	if(set_u_init)
		{
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_out_set(config, dims, out, ii, "u", u_init+ii*nu);
			}
		}
	else // initialize to zero
		{
		double *u_init = calloc(nu, sizeof(double));
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_out_set(config, dims, out, ii, "u", u_init);
			}
		free(u_init);
		}



	/* solver */

	ocp_nlp_solver *solver = ocp_nlp_solver_create(config, dims, opts);



	/* populate output struct */

	// config
	mxArray *config_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
	l_ptr = mxGetData(config_mat);
	l_ptr[0] = (long long) config;
	mxSetField(plhs[0], 0, "config", config_mat);

	// dims
	mxArray *dims_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
	l_ptr = mxGetData(dims_mat);
	l_ptr[0] = (long long) dims;
	mxSetField(plhs[0], 0, "dims", dims_mat);

	// opts
	mxArray *opts_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
	l_ptr = mxGetData(opts_mat);
	l_ptr[0] = (long long) opts;
	mxSetField(plhs[0], 0, "opts", opts_mat);

	// in
	mxArray *in_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
	l_ptr = mxGetData(in_mat);
	l_ptr[0] = (long long) in;
	mxSetField(plhs[0], 0, "in", in_mat);

	// out
	mxArray *out_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
	l_ptr = mxGetData(out_mat);
	l_ptr[0] = (long long) out;
	mxSetField(plhs[0], 0, "out", out_mat);

	// solver
	mxArray *solver_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
	l_ptr = mxGetData(solver_mat);
	l_ptr[0] = (long long) solver;
	mxSetField(plhs[0], 0, "solver", solver_mat);



	/* return */
	return;

	}

