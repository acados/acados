// system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
// acados
#include "acados_c/ocp_nlp_interface.h"
#include "acados/utils/external_function_generic.h"
#include "acados_c/external_function_interface.h"
// mex
#include "mex.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
	{

//	mexPrintf("\nin ocp_ext_fun_destroy\n");

	int ii;
	long long *ptr;



	/* RHS */

	// model

	// TODO bool instead !!!
	char *param_f = mxArrayToString( mxGetField( prhs[0], 0, "dyn_param_f" ) );
	char *param_phi = mxArrayToString( mxGetField( prhs[0], 0, "dyn_param_phi" ) );
	char *param_h = mxArrayToString( mxGetField( prhs[0], 0, "constr_param_h" ) );
	char *param_h_e = mxArrayToString( mxGetField( prhs[0], 0, "constr_param_h_e" ) );
	char *param_y = mxArrayToString( mxGetField( prhs[0], 0, "cost_param_y" ) );
	char *param_y_e = mxArrayToString( mxGetField( prhs[0], 0, "cost_param_y_e" ) );
	char *param_ext_cost = mxArrayToString( mxGetField( prhs[0], 0, "cost_param_ext_cost" ) );
	char *param_ext_cost_e = mxArrayToString( mxGetField( prhs[0], 0, "cost_param_ext_cost_e" ) );

	// C_ocp

	// dims
	ptr = (long long *) mxGetData( mxGetField( prhs[1], 0, "dims" ) );
	ocp_nlp_dims *dims = (ocp_nlp_dims *) ptr[0];

	int N = dims->N;


	// C_ocp_ext_fun

	external_function_casadi *ext_fun_ptr;
	external_function_param_casadi *ext_fun_param_ptr;

	// TODO use a loop and mxGetFeildByNumber ?
	// NO! it needs to know if a fun refers to middle or last stage !!!

	if(mxGetField( prhs[2], 0, "dyn_expl_ode_fun" )!=NULL && mxGetM(mxGetField( prhs[2], 0, "dyn_expl_ode_fun" ))>0)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "dyn_expl_ode_fun" ) );
		if(!strcmp(param_f, "true")) // TODO bool
			{
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				external_function_param_casadi_free(ext_fun_param_ptr+ii);
				}
			free(ext_fun_param_ptr);
			}
		else
			{
			ext_fun_ptr = (external_function_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				external_function_casadi_free(ext_fun_ptr+ii);
				}
			free(ext_fun_ptr);
			}
		}
	if(mxGetField( prhs[2], 0, "dyn_expl_vde_for" )!=NULL && mxGetM(mxGetField( prhs[2], 0, "dyn_expl_vde_for" ))>0)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "dyn_expl_vde_for" ) );
		if(!strcmp(param_f, "true")) // TODO bool
			{
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				external_function_param_casadi_free(ext_fun_param_ptr+ii);
				}
			free(ext_fun_param_ptr);
			}
		else
			{
			ext_fun_ptr = (external_function_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				external_function_casadi_free(ext_fun_ptr+ii);
				}
			free(ext_fun_ptr);
			}
		}
	if(mxGetField( prhs[2], 0, "dyn_expl_vde_adj" )!=NULL && mxGetM(mxGetField( prhs[2], 0, "dyn_expl_vde_adj" ))>0)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "dyn_expl_vde_adj" ) );
		if(!strcmp(param_f, "true")) // TODO bool
			{
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				external_function_param_casadi_free(ext_fun_param_ptr+ii);
				}
			free(ext_fun_param_ptr);
			}
		else
			{
			ext_fun_ptr = (external_function_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				external_function_casadi_free(ext_fun_ptr+ii);
				}
			free(ext_fun_ptr);
			}
		}
	if(mxGetField( prhs[2], 0, "dyn_expl_ode_hes" )!=NULL && mxGetM(mxGetField( prhs[2], 0, "dyn_expl_ode_hes" ))>0)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "dyn_expl_ode_hes" ) );
		if(!strcmp(param_f, "true")) // TODO bool
			{
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				external_function_param_casadi_free(ext_fun_param_ptr+ii);
				}
			free(ext_fun_param_ptr);
			}
		else
			{
			ext_fun_ptr = (external_function_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				external_function_casadi_free(ext_fun_ptr+ii);
				}
			free(ext_fun_ptr);
			}
		}
	if(mxGetField( prhs[2], 0, "dyn_impl_ode_fun" )!=NULL && mxGetM(mxGetField( prhs[2], 0, "dyn_impl_ode_fun" ))>0)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "dyn_impl_ode_fun" ) );
		if(!strcmp(param_f, "true")) // TODO bool
			{
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				external_function_param_casadi_free(ext_fun_param_ptr+ii);
				}
			free(ext_fun_param_ptr);
			}
		else
			{
			ext_fun_ptr = (external_function_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				external_function_casadi_free(ext_fun_ptr+ii);
				}
			free(ext_fun_ptr);
			}
		}
	if(mxGetField( prhs[2], 0, "dyn_impl_ode_fun_jac_x_xdot" )!=NULL && mxGetM(mxGetField( prhs[2], 0, "dyn_impl_ode_fun_jac_x_xdot" ))>0)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "dyn_impl_ode_fun_jac_x_xdot" ) );
		if(!strcmp(param_f, "true")) // TODO bool
			{
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				external_function_param_casadi_free(ext_fun_param_ptr+ii);
				}
			free(ext_fun_param_ptr);
			}
		else
			{
			ext_fun_ptr = (external_function_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				external_function_casadi_free(ext_fun_ptr+ii);
				}
			free(ext_fun_ptr);
			}
		}
	if(mxGetField( prhs[2], 0, "dyn_impl_ode_jac_x_xdot_u" )!=NULL && mxGetM(mxGetField( prhs[2], 0, "dyn_impl_ode_jac_x_xdot_u" ))>0)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "dyn_impl_ode_jac_x_xdot_u" ) );
		if(!strcmp(param_f, "true")) // TODO bool
			{
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				external_function_param_casadi_free(ext_fun_param_ptr+ii);
				}
			free(ext_fun_param_ptr);
			}
		else
			{
			ext_fun_ptr = (external_function_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				external_function_casadi_free(ext_fun_ptr+ii);
				}
			free(ext_fun_ptr);
			}
		}
	if(mxGetField( prhs[2], 0, "dyn_impl_ode_hess" )!=NULL && mxGetM(mxGetField( prhs[2], 0, "dyn_impl_ode_hess" ))>0)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "dyn_impl_ode_hess" ) );
		if(!strcmp(param_f, "true")) // TODO bool
			{
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				external_function_param_casadi_free(ext_fun_param_ptr+ii);
				}
			free(ext_fun_param_ptr);
			}
		else
			{
			ext_fun_ptr = (external_function_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				external_function_casadi_free(ext_fun_ptr+ii);
				}
			free(ext_fun_ptr);
			}
		}
	if(mxGetField( prhs[2], 0, "dyn_disc_phi_fun_jac" )!=NULL && mxGetM(mxGetField( prhs[2], 0, "dyn_disc_phi_fun_jac" ))>0)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "dyn_disc_phi_fun_jac" ) );
		if(!strcmp(param_f, "true")) // TODO bool
			{
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				external_function_param_casadi_free(ext_fun_param_ptr+ii);
				}
			free(ext_fun_param_ptr);
			}
		else
			{
			ext_fun_ptr = (external_function_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				external_function_casadi_free(ext_fun_ptr+ii);
				}
			free(ext_fun_ptr);
			}
		}
	if(mxGetField( prhs[2], 0, "dyn_disc_phi_fun_jac_hess" )!=NULL && mxGetM(mxGetField( prhs[2], 0, "dyn_disc_phi_fun_jac_hess" ))>0)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "dyn_disc_phi_fun_jac_hess" ) );
		if(!strcmp(param_f, "true")) // TODO bool
			{
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				external_function_param_casadi_free(ext_fun_param_ptr+ii);
				}
			free(ext_fun_param_ptr);
			}
		else
			{
			ext_fun_ptr = (external_function_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				external_function_casadi_free(ext_fun_ptr+ii);
				}
			free(ext_fun_ptr);
			}
		}
	if(mxGetField( prhs[2], 0, "constr_h_fun_jac_ut_xt" )!=NULL && mxGetM(mxGetField( prhs[2], 0, "constr_h_fun_jac_ut_xt" ))>0)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "constr_h_fun_jac_ut_xt" ) );
		if(!strcmp(param_h, "true")) // TODO bool
			{
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				external_function_param_casadi_free(ext_fun_param_ptr+ii);
				}
			free(ext_fun_param_ptr);
			}
		else
			{
			ext_fun_ptr = (external_function_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				external_function_casadi_free(ext_fun_ptr+ii);
				}
			free(ext_fun_ptr);
			}
		}
	if(mxGetField( prhs[2], 0, "constr_h_fun_jac_ut_xt_hess" )!=NULL && mxGetM(mxGetField( prhs[2], 0, "constr_h_fun_jac_ut_xt_hess" ))>0)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "constr_h_fun_jac_ut_xt_hess" ) );
		if(!strcmp(param_h, "true")) // TODO bool
			{
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				external_function_param_casadi_free(ext_fun_param_ptr+ii);
				}
			free(ext_fun_param_ptr);
			}
		else
			{
			ext_fun_ptr = (external_function_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				external_function_casadi_free(ext_fun_ptr+ii);
				}
			free(ext_fun_ptr);
			}
		}
	if(mxGetField( prhs[2], 0, "constr_h_e_fun_jac_ut_xt" )!=NULL)
	if(mxGetField( prhs[2], 0, "constr_h_e_fun_jac_ut_xt" )!=NULL && mxGetM(mxGetField( prhs[2], 0, "constr_h_e_fun_jac_ut_xt" ))>0)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "constr_h_e_fun_jac_ut_xt" ) );
		if(!strcmp(param_h_e, "true")) // TODO bool
			{
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			external_function_param_casadi_free(ext_fun_param_ptr);
			free(ext_fun_param_ptr);
			}
		else
			{
			ext_fun_ptr = (external_function_casadi *) ptr[0];
			external_function_casadi_free(ext_fun_ptr);
			free(ext_fun_ptr);
			}
		}
	if(mxGetField( prhs[2], 0, "constr_h_e_fun_jac_ut_xt_hess" )!=NULL && mxGetM(mxGetField( prhs[2], 0, "constr_h_e_fun_jac_ut_xt_hess" ))>0)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "constr_h_e_fun_jac_ut_xt_hess" ) );
		if(!strcmp(param_h_e, "true")) // TODO bool
			{
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			external_function_param_casadi_free(ext_fun_param_ptr);
			free(ext_fun_param_ptr);
			}
		else
			{
			ext_fun_ptr = (external_function_casadi *) ptr[0];
			external_function_casadi_free(ext_fun_ptr);
			free(ext_fun_ptr);
			}
		}
	if(mxGetField( prhs[2], 0, "cost_y_fun_jac_ut_xt" )!=NULL && mxGetM(mxGetField( prhs[2], 0, "cost_y_fun_jac_ut_xt" ))>0)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "cost_y_fun_jac_ut_xt" ) );
		if(!strcmp(param_y, "true")) // TODO bool
			{
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				external_function_param_casadi_free(ext_fun_param_ptr+ii);
				}
			free(ext_fun_param_ptr);
			}
		else
			{
			ext_fun_ptr = (external_function_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				external_function_casadi_free(ext_fun_ptr+ii);
				}
			free(ext_fun_ptr);
			}
		}
	if(mxGetField( prhs[2], 0, "cost_y_hess" )!=NULL && mxGetM(mxGetField( prhs[2], 0, "cost_y_hess" ))>0)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "cost_y_hess" ) );
		if(!strcmp(param_y, "true")) // TODO bool
			{
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				external_function_param_casadi_free(ext_fun_param_ptr+ii);
				}
			free(ext_fun_param_ptr);
			}
		else
			{
			ext_fun_ptr = (external_function_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				external_function_casadi_free(ext_fun_ptr+ii);
				}
			free(ext_fun_ptr);
			}
		}
	if(mxGetField( prhs[2], 0, "cost_y_e_fun_jac_ut_xt" )!=NULL && mxGetM(mxGetField( prhs[2], 0, "cost_y_e_fun_jac_ut_xt" ))>0)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "cost_y_e_fun_jac_ut_xt" ) );
		if(!strcmp(param_y_e, "true")) // TODO bool
			{
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			external_function_param_casadi_free(ext_fun_param_ptr);
			free(ext_fun_param_ptr);
			}
		else
			{
			ext_fun_ptr = (external_function_casadi *) ptr[0];
			external_function_casadi_free(ext_fun_ptr);
			free(ext_fun_ptr);
			}
		}
	if(mxGetField( prhs[2], 0, "cost_y_e_hess" )!=NULL && mxGetM(mxGetField( prhs[2], 0, "cost_y_e_hess" ))>0)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "cost_y_e_hess" ) );
		if(!strcmp(param_y_e, "true")) // TODO bool
			{
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			external_function_param_casadi_free(ext_fun_param_ptr);
			free(ext_fun_param_ptr);
			}
		else
			{
			ext_fun_ptr = (external_function_casadi *) ptr[0];
			external_function_casadi_free(ext_fun_ptr);
			free(ext_fun_ptr);
			}
		}
	if(mxGetField( prhs[2], 0, "cost_ext_cost_jac_hes" )!=NULL && mxGetM(mxGetField( prhs[2], 0, "cost_ext_cost_jac_hes" ))>0)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "cost_ext_cost_jac_hes" ) );
		if(!strcmp(param_ext_cost, "true")) // TODO bool
			{
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				external_function_param_casadi_free(ext_fun_param_ptr+ii);
				}
			free(ext_fun_param_ptr);
			}
		else
			{
			ext_fun_ptr = (external_function_casadi *) ptr[0];
			for(ii=0; ii<N; ii++)
				{
				external_function_casadi_free(ext_fun_ptr+ii);
				}
			free(ext_fun_ptr);
			}
		}
	if(mxGetField( prhs[2], 0, "cost_ext_cost_e_jac_hes" )!=NULL && mxGetM(mxGetField( prhs[2], 0, "cost_ext_cost_e_jac_hes" ))>0)
		{
		ptr = (long long *) mxGetData( mxGetField( prhs[2], 0, "cost_ext_cost_e_jac_hes" ) );
		if(!strcmp(param_ext_cost_e, "true")) // TODO bool
			{
			ext_fun_param_ptr = (external_function_param_casadi *) ptr[0];
			external_function_param_casadi_free(ext_fun_param_ptr);
			free(ext_fun_param_ptr);
			}
		else
			{
			ext_fun_ptr = (external_function_casadi *) ptr[0];
			external_function_casadi_free(ext_fun_ptr);
			free(ext_fun_ptr);
			}
		}



	/* return */

	return;

	}


