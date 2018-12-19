%% create external cost function for nonlinear chain model
% import casadi
% this script was written to replace hand written external functions that
% looked as pasted in the end of the script
clearvars

import casadi.*
addpath('../../../interfaces/matlab/external_function_generation');

for nm = 2:6
    source = strcat('xN_nm', num2str(nm), '.txt');
    xref = load(source);
    nx = 6 * (nm-1);
    nu = 3;
    uxref = zeros(nu + nx, 1);
    for i = 1:nx
        uxref(i+nu) = xref(i);
    end

    model.x = SX.sym('x',nx,1);
    model.u = SX.sym('u',nu,1);
    cost.name = strcat('chain_nm_', num2str(nm));
    
    ux = [model.u; model.x];

    % define least square cost as general expression
    cost.general_expr = SX.zeros(1,1);
    for j = 1:nu
        cost.general_expr = cost.general_expr + 0.5 * (ux(j) - uxref(j))^2;
    end
    for j = 1:nx
        cost.general_expr = cost.general_expr + 0.5 * 1e-2 * (ux(j+nu) - uxref(j+nu))^2;
    end
    
    generate_c_code_external_cost(model, cost)
    
end


% 
% // hand-generated external function for externally provided hessian and gradient
% void ext_cost_nm2(void *fun, ext_fun_arg_t *type_in, void **in, ext_fun_arg_t *type_out, void **out)
% {
% 
% 	int ii;
% 
% 	int nu = 3;
% 	int nx = 6;
% 
% 	int nv = nu+nx;
% 
% 	// ref
% 	double *ref = calloc(nx+nu, sizeof(double));
% 	for (ii=0; ii<nu; ii++)
% 		ref[ii] = 0.0;
% 	for (ii=0; ii<nx; ii++)
% 		ref[nu+ii] = xN_nm2[ii];
% 
% 	// Hessian
% 	double *hess = out[1];
% 	for (ii=0; ii<nv*nv; ii++)
% 		hess[ii] = 0.0;
% 	for (ii=0; ii<nu; ii++)
% 		hess[ii*(nv+1)] = 1.0;
% 	for (; ii<nu+nx; ii++)
% 		hess[ii*(nv+1)] = 1e-2;
% 
% 	// gradient
% 	double *ux = in[0];
% 	double *grad = out[0];
% 	for (ii=0; ii<nv; ii++)
% 		grad[ii] = 0.0;
% 	for (ii=0; ii<nv; ii++)
% 		grad[ii] = hess[ii*(nv+1)] * (ux[ii] - ref[ii]);
% 
%     free(ref);
% 	return;
% 
% }
% 
% void ext_cost_nm3(void *fun, ext_fun_arg_t *type_in, void **in, ext_fun_arg_t *type_out, void **out)
% {
% 
% 	int ii;
% 
% 	int nu = 3;
% 	int nx = 12;
% 
% 	int nv = nu+nx;
% 
% 	// ref
%     double *ref = calloc(nx+nu, sizeof(double));
% 	for (ii=0; ii<nu; ii++)
% 		ref[ii] = 0.0;
% 	for (ii=0; ii<nx; ii++)
% 		ref[nu+ii] = xN_nm3[ii];
% 
% 	// Hessian
% 	double *hess = out[1];
% 	for (ii=0; ii<nv*nv; ii++)
% 		hess[ii] = 0.0;
% 	for (ii=0; ii<nu; ii++)
% 		hess[ii*(nv+1)] = 1.0;
% 	for (; ii<nu+nx; ii++)
% 		hess[ii*(nv+1)] = 1e-2;
% 
% 	// gradient
% 	double *ux = in[0];
% 	double *grad = out[0];
% 	for (ii=0; ii<nv; ii++)
% 		grad[ii] = 0.0;
% 	for (ii=0; ii<nv; ii++)
% 		grad[ii] = hess[ii*(nv+1)] * (ux[ii] - ref[ii]);
% 
%     free(ref);
% 	return;
% 
% }
