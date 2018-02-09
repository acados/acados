// standard
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
// acados
#include "acados/utils/print.h"
#include "acados/utils/mem.h"

#include "acados/sim/sim_common.h"
#include "acados/sim/sim_irk_integrator.h"

#include "external/blasfeo/include/blasfeo_target.h"
#include "external/blasfeo/include/blasfeo_common.h"
#include "external/blasfeo/include/blasfeo_d_aux.h"
#include "external/blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "external/blasfeo/include/blasfeo_v_aux_ext_dep.h"
#include "external/blasfeo/include/blasfeo_d_blas.h"



#define IMPL_RES_NUMIN 4
#define IMPL_RES_NUMOUT 1
#define IMPL_JAC_NUMIN 4
#define IMPL_JAC_NUMOUT 3



int sim_irk_impl_res_input_dims[IMPL_RES_NUMIN] = {0};
int sim_irk_impl_res_output_dims[IMPL_RES_NUMOUT] = {0};
external_function_dims sim_irk_impl_res_dims = {
    IMPL_RES_NUMIN,
    IMPL_RES_NUMOUT,
    sim_irk_impl_res_input_dims,
    sim_irk_impl_res_output_dims
};



int sim_irk_impl_jac_input_dims[IMPL_JAC_NUMIN] = {0};
int sim_irk_impl_jac_output_dims[IMPL_JAC_NUMOUT] = {0};
external_function_dims sim_irk_impl_jac_dims = {
    IMPL_JAC_NUMIN,
    IMPL_JAC_NUMOUT,
    sim_irk_impl_jac_input_dims,
    sim_irk_impl_jac_output_dims
};



int sim_irk_integrator_calculate_args_size(sim_dims *dims, void *submodules_)
{
    sim_irk_integrator_submodules *submodules = (sim_irk_integrator_submodules *) submodules_;

    int size = sizeof(sim_irk_integrator_args);

    int ns = dims->num_stages;
    size += ns * ns * sizeof(double);  // A_mat
    size += ns * sizeof(double);  // b_vec
    size += ns * sizeof(double);  // c_vec

    size += 2*sizeof(external_function_fcn_ptrs);

    if (submodules->impl_res != NULL) {
        size += submodules->impl_res->calculate_args_size(&sim_irk_impl_res_dims, submodules->impl_res->submodules);
    }
    
    if (submodules->impl_jac != NULL) {
        size += submodules->impl_jac->calculate_args_size(&sim_irk_impl_jac_dims, submodules->impl_jac->submodules);
    }

    make_int_multiple_of(8, &size);
    size += 1 * 8;

    return size;
}



void *sim_irk_integrator_assign_args(sim_dims *dims, void **submodules_, void *raw_memory)
{
    sim_irk_integrator_submodules *submodules = (sim_irk_integrator_submodules *) *submodules_;

    char *c_ptr = (char *) raw_memory;

    sim_irk_integrator_args *args = (sim_irk_integrator_args *) c_ptr;
    c_ptr += sizeof(sim_irk_integrator_args);

    int ns = dims->num_stages;
    args->num_stages = ns;

    align_char_to(8, &c_ptr);

    assign_double(ns*ns, &args->A_mat, &c_ptr);
    assign_double(ns, &args->b_vec, &c_ptr);
    assign_double(ns, &args->c_vec, &c_ptr);

    if (submodules->impl_res != NULL) {
        args->submodules.impl_res = (external_function_fcn_ptrs *)c_ptr;
        c_ptr += sizeof(external_function_fcn_ptrs);

        void *impl_res_submodules = submodules->impl_res->submodules;
        args->impl_res_args = submodules->impl_res->assign_args(&sim_irk_impl_res_dims, &(impl_res_submodules), c_ptr);
        c_ptr += submodules->impl_res->calculate_args_size(&sim_irk_impl_res_dims, submodules->impl_res->submodules);

        *(args->submodules.impl_res) = *(submodules->impl_res);
        args->submodules.impl_res->submodules = impl_res_submodules;
    } else {
        args->submodules.impl_res = NULL;
    }

    if (submodules->impl_jac != NULL) {
        args->submodules.impl_jac = (external_function_fcn_ptrs *)c_ptr;
        c_ptr += sizeof(external_function_fcn_ptrs);

        void *impl_jac_submodules = submodules->impl_jac->submodules;
        args->impl_jac_args = submodules->impl_jac->assign_args(&sim_irk_impl_jac_dims, &(impl_jac_submodules), c_ptr);
        c_ptr += submodules->impl_jac->calculate_args_size(&sim_irk_impl_jac_dims, submodules->impl_jac->submodules);

        *(args->submodules.impl_jac) = *(submodules->impl_jac);
        args->submodules.impl_jac->submodules = impl_jac_submodules;
    } else {
        args->submodules.impl_jac = NULL;
    }

    assert((char*)raw_memory + sim_irk_integrator_calculate_args_size(dims, *submodules_) >= c_ptr);

    args->newton_iter = 3;

    // Update submodules pointer
    *submodules_ = (void *)&(args->submodules);

    return (void *)args;
}



void *sim_irk_integrator_copy_args(sim_dims *dims, void *raw_memory, void *source_)
{
    sim_irk_integrator_args *source = (sim_irk_integrator_args *)source_;
    sim_irk_integrator_args *dest;

    sim_irk_integrator_submodules *submodules = &source->submodules;

    dest = sim_irk_integrator_assign_args(dims, (void **) &submodules, raw_memory);

    dest->interval = source->interval;
    dest->num_stages = source->num_stages;
    dest->num_steps = source->num_steps;
    dest->num_forw_sens = source->num_forw_sens;
    dest->sens_forw = source->sens_forw;
    dest->sens_adj = source->sens_adj;
    dest->sens_hess = source->sens_hess;
    dest->newton_iter = source->newton_iter;
    dest->jac_reuse = source->jac_reuse;

    int ns = dims->num_stages;

    memcpy(dest->A_mat, source->A_mat, ns*ns*sizeof(double));
    memcpy(dest->c_vec, source->c_vec, ns*sizeof(double));
    memcpy(dest->b_vec, source->b_vec, ns*sizeof(double));

    source->submodules.impl_res->copy_args(&sim_irk_impl_res_dims, dest->impl_res_args, source->impl_res_args);
    source->submodules.impl_jac->copy_args(&sim_irk_impl_jac_dims, dest->impl_jac_args, source->impl_jac_args);

    return (void *)dest;
}



void sim_irk_integrator_initialize_default_args(sim_dims *dims, void *args_)
{
    sim_irk_integrator_args *args = (sim_irk_integrator_args *) args_;
    int ns = args->num_stages;

    assert(args->num_stages == 3 && "only number of stages = 3 implemented!");

    memcpy(args->A_mat,
        ((real_t[]){0.1389, 0.3003, 0.2680 ,
        -0.0360, 0.2222, 0.4804,
        0.0098, -0.0225, 0.1389}),
        sizeof(*args->A_mat) * (ns * ns));
    memcpy(args->b_vec, ((real_t[]){0.2778, 0.4444, 0.2778}),
        sizeof(*args->b_vec) * (ns));
    memcpy(args->c_vec, ((real_t[]){0.1127, 0.5, 0.8873}),
        sizeof(*args->c_vec) * (ns));

    args->num_steps = 2;
    args->num_forw_sens = dims->nx + dims->nu;
    args->sens_forw = true;
    args->sens_adj = false;
    args->sens_hess = false;
    args->jac_reuse = false;

    args->submodules.impl_res->initialize_default_args(args->impl_res_args);
    args->submodules.impl_jac->initialize_default_args(args->impl_jac_args);
}



int sim_irk_integrator_calculate_memory_size(sim_dims *dims, void *args_)
{
    sim_irk_integrator_args *args = (sim_irk_integrator_args *)args_;

    int size = sizeof(sim_irk_integrator_memory);

    size += args->submodules.impl_res->calculate_memory_size(&sim_irk_impl_res_dims, args->impl_res_args);
    
    size += args->submodules.impl_jac->calculate_memory_size(&sim_irk_impl_jac_dims, args->impl_jac_args);

    make_int_multiple_of(8, &size);
    size += 1 * 8;

    return size;
}



void *sim_irk_integrator_assign_memory(sim_dims *dims, void *args_, void *raw_memory)
{
    sim_irk_integrator_args *args = (sim_irk_integrator_args *)args_;

    sim_irk_integrator_memory *mem;

    char *c_ptr = (char *) raw_memory;

    mem = (sim_irk_integrator_memory *) c_ptr;
    c_ptr += sizeof(sim_irk_integrator_memory);

    mem->impl_res_mem = args->submodules.impl_res->assign_memory(&sim_irk_impl_res_dims, args->impl_res_args, c_ptr);
    c_ptr += args->submodules.impl_res->calculate_memory_size(&sim_irk_impl_res_dims, args->impl_res_args);

    mem->impl_jac_mem = args->submodules.impl_jac->assign_memory(&sim_irk_impl_jac_dims, args->impl_jac_args, c_ptr);
    c_ptr += args->submodules.impl_jac->calculate_memory_size(&sim_irk_impl_jac_dims, args->impl_jac_args);

    align_char_to(8, &c_ptr);

    assert((size_t)c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    assert((char*)raw_memory + sim_irk_integrator_calculate_memory_size(dims, args_) >= c_ptr);

    return (void *)mem;
}



int sim_irk_integrator_calculate_workspace_size(sim_dims *dims, void *args_)
{
    sim_irk_integrator_args *args = (sim_irk_integrator_args *) args_;

    int nx = dims->nx;
    int nu = dims->nu;
    int np = dims->np;
    
    int ns = args->num_stages; // number of stages
    int steps = args->num_steps;

    int size = sizeof(sim_irk_integrator_workspace);

    size += 4*sizeof(struct blasfeo_dmat); // JGK, JGf, JKf, S_forw
    size += steps*sizeof(struct blasfeo_dmat); // JG_traj

    size += 4*sizeof(struct blasfeo_dvec); // rG, K, xt, xn
    size += 2*sizeof(struct blasfeo_dvec); // lambda,lambdaK
    size += 2*steps*sizeof(struct blasfeo_dvec);  // **xn_traj, **K_traj;
    
    size += blasfeo_memsize_dmat(nx*ns, nx*ns); // JGK
    size += 2*blasfeo_memsize_dmat(nx*ns, nx+nu); // JGf, JKf
    size += blasfeo_memsize_dmat(nx, nx+nu); // S_forw
    size += steps * blasfeo_memsize_dmat(nx*ns, nx*ns); // for JG_traj  

    size += 2*blasfeo_memsize_dvec(nx*ns); // rG, K
    size += 2*blasfeo_memsize_dvec(nx); // xt, x
    size += blasfeo_memsize_dvec(nx+nu); // lambda
    size += blasfeo_memsize_dvec(nx*ns); // lambdaK
    size += steps * blasfeo_memsize_dvec(nx); // for xn_traj
    size += steps * blasfeo_memsize_dvec(nx*ns); // for K_traj
    
    size += nx * sizeof(double); //  rGt
    size += nx * (2*nx+nu) * sizeof(double); // jac_out
    size += nx * nx * sizeof(double); // Jt
    size += (2*nx + nu + np) * sizeof(double); // ode_args
    size += (nx+nu) * sizeof(double); // S_adj_w

    size += nx *ns * (steps+1) * sizeof(int); // ipiv

    size += args->submodules.impl_res->calculate_workspace_size(&sim_irk_impl_res_dims, args->impl_res_args);

    size += args->submodules.impl_jac->calculate_workspace_size(&sim_irk_impl_jac_dims, args->impl_jac_args);
    
    make_int_multiple_of(64, &size);
    size += 1 * 64; 

    return size;
}



static void *cast_workspace(sim_dims *dims, void *args_, void *raw_memory)
{
    sim_irk_integrator_args *args = (sim_irk_integrator_args *) args_;

    int nx = dims->nx;
    int nu = dims->nu;
    int np = dims->np;
        
    int ns = args->num_stages;
    int steps = args->num_steps;

    char *c_ptr = (char *)raw_memory;

    sim_irk_integrator_workspace *workspace = (sim_irk_integrator_workspace *) c_ptr;
    c_ptr += sizeof(sim_irk_integrator_workspace);

    assign_blasfeo_dmat_structs(steps, &workspace->JG_traj, &c_ptr);

    workspace->JGK = (struct blasfeo_dmat *)c_ptr;
    c_ptr += sizeof(struct blasfeo_dmat); 

    workspace->JGf = (struct blasfeo_dmat *)c_ptr;
    c_ptr += sizeof(struct blasfeo_dmat);

    workspace->JKf = (struct blasfeo_dmat *)c_ptr;
    c_ptr += sizeof(struct blasfeo_dmat);

    workspace->S_forw = (struct blasfeo_dmat *)c_ptr;
    c_ptr += sizeof(struct blasfeo_dmat);


    assign_blasfeo_dvec_structs(steps, &workspace->xn_traj, &c_ptr);
    assign_blasfeo_dvec_structs(steps, &workspace->K_traj, &c_ptr);

    workspace->rG = (struct blasfeo_dvec *)c_ptr;
    c_ptr += sizeof(struct blasfeo_dvec);

    workspace->K = (struct blasfeo_dvec *)c_ptr;
    c_ptr += sizeof(struct blasfeo_dvec);

    workspace->xt = (struct blasfeo_dvec *)c_ptr;
    c_ptr += sizeof(struct blasfeo_dvec);

    workspace->xn = (struct blasfeo_dvec *)c_ptr;
    c_ptr += sizeof(struct blasfeo_dvec);

    workspace->lambda = (struct blasfeo_dvec *)c_ptr;
    c_ptr += sizeof(struct blasfeo_dvec);

    workspace->lambdaK = (struct blasfeo_dvec *)c_ptr;
    c_ptr += sizeof(struct blasfeo_dvec);
    
    align_char_to(64, &c_ptr);

    assign_blasfeo_dmat_mem(nx*ns, nx*ns, workspace->JGK, &c_ptr);
    assign_blasfeo_dmat_mem(nx*ns, nx+nu, workspace->JGf, &c_ptr);
    assign_blasfeo_dmat_mem(nx*ns, nx+nu, workspace->JKf, &c_ptr);
    assign_blasfeo_dmat_mem(nx, nx+nu, workspace->S_forw, &c_ptr);
    for (int i=0;i<steps;i++){
        assign_blasfeo_dmat_mem(nx*ns, nx*ns, &workspace->JG_traj[i], &c_ptr);
    }

    assign_blasfeo_dvec_mem(nx*ns, workspace->rG, &c_ptr);
    assign_blasfeo_dvec_mem(nx*ns, workspace->K, &c_ptr);
    assign_blasfeo_dvec_mem(nx, workspace->xt, &c_ptr);
    assign_blasfeo_dvec_mem(nx, workspace->xn, &c_ptr);
    assign_blasfeo_dvec_mem(nx+nu, workspace->lambda, &c_ptr);
    assign_blasfeo_dvec_mem(nx*ns, workspace->lambdaK, &c_ptr);
    for (int i=0;i<steps;i++){
        assign_blasfeo_dvec_mem(nx, &workspace->xn_traj[i], &c_ptr);
        assign_blasfeo_dvec_mem(nx*ns, &workspace->K_traj[i], &c_ptr);    
    }

    assign_double(nx, &workspace->rGt, &c_ptr);
    assign_double(nx * (2*nx+nu), &workspace->jac_out, &c_ptr);
    assign_double(nx * nx, &workspace->Jt, &c_ptr);
    assign_double(2*nx + nu + np, &workspace->ode_args, &c_ptr);
    assign_double(nx + nu, &workspace->S_adj_w, &c_ptr);

    assign_int(nx * ns* (steps+1) , &workspace->ipiv, &c_ptr);

    workspace->impl_res_work = (void *) c_ptr;
    c_ptr += args->submodules.impl_res->calculate_workspace_size(&sim_irk_impl_res_dims, args->impl_res_args);

    workspace->impl_jac_work = (void *) c_ptr;
    c_ptr += args->submodules.impl_jac->calculate_workspace_size(&sim_irk_impl_jac_dims, args->impl_jac_args);

    // printf("\npointer moved - size calculated = %d bytes\n", c_ptr- (char*)raw_memory - sim_irk_integrator_calculate_workspace_size(dims, args_));

    assert((char*)raw_memory + sim_irk_integrator_calculate_workspace_size(dims, args_) >= c_ptr);

    return (void *)workspace;
}



static void compute_impl_res(const int nx, const int nu, const int np, double *in, double *out, 
                         sim_irk_integrator_args *args,
                         sim_irk_integrator_memory *mem,
                         sim_irk_integrator_workspace *work) 
{
    double *x = in;
    double *xdot = in + nx;
    double *u = in + 2 * nx;
    double *p = in + 2 * nx + nu;

    double *res = out;

    double *impl_res_in_inputs[IMPL_RES_NUMIN];
    impl_res_in_inputs[0] = x;
    impl_res_in_inputs[1] = xdot;
    impl_res_in_inputs[2] = u;
    impl_res_in_inputs[3] = p;

    bool impl_res_in_compute_output[IMPL_RES_NUMOUT] = {true};

    double *impl_res_out_outputs[IMPL_RES_NUMOUT];
    impl_res_out_outputs[0] = res;

    external_function_in impl_res_in;
    impl_res_in.inputs = impl_res_in_inputs;
    impl_res_in.compute_output = impl_res_in_compute_output;

    external_function_out impl_res_out;
    impl_res_out.outputs = impl_res_out_outputs;

    args->submodules.impl_res->fun(&impl_res_in, &impl_res_out, args->impl_res_args, mem->impl_res_mem, work->impl_res_work);
}



static void compute_impl_jac(const int nx, const int nu, const int np, double *in, double *out, 
                         sim_irk_integrator_args *args,
                         sim_irk_integrator_memory *mem,
                         sim_irk_integrator_workspace *work)
{
    double *x = in;
    double *xdot = in + nx;
    double *u = in + 2 * nx;
    double *p = in + 2 * nx + nu;

    double *jac_x = out;
    double *jac_xdot = out + nx*nx;
    double *jac_u = out + 2*nx*nx;

    double *impl_jac_in_inputs[IMPL_JAC_NUMIN];
    impl_jac_in_inputs[0] = x;
    impl_jac_in_inputs[1] = xdot;
    impl_jac_in_inputs[2] = u;
    impl_jac_in_inputs[3] = p;

    bool impl_jac_in_compute_output[IMPL_JAC_NUMOUT] = {true, true, true};

    double *impl_jac_out_outputs[IMPL_JAC_NUMOUT];
    impl_jac_out_outputs[0] = jac_x;
    impl_jac_out_outputs[1] = jac_xdot;
    impl_jac_out_outputs[2] = jac_u;

    external_function_in impl_jac_in;
    impl_jac_in.inputs = impl_jac_in_inputs;
    impl_jac_in.compute_output = impl_jac_in_compute_output;

    external_function_out impl_jac_out;
    impl_jac_out.outputs = impl_jac_out_outputs;

    args->submodules.impl_jac->fun(&impl_jac_in, &impl_jac_out, args->impl_jac_args, mem->impl_jac_mem, work->impl_jac_work);
}



int sim_irk_integrator(sim_in *in, sim_out *out, void *args_, void *mem_, void *work_)
{
    sim_irk_integrator_args *args = (sim_irk_integrator_args *) args_;
    sim_irk_integrator_memory *mem = (sim_irk_integrator_memory *) mem_;

    sim_dims dims = {
        args->num_stages,
        in->nx,
        in->nu,
        in->np
    };
    sim_irk_integrator_workspace *workspace = (sim_irk_integrator_workspace *) cast_workspace(&dims, args, work_);

    int ii, jj, iter, kk, ss;
    double a,b;

    int nx = in->nx;
	int nu = in->nu;
	int np = in->np;
    int num_steps = args->num_steps;
    double *x = in->x;
	double *u = in->u;
    double *p = in->p;
    double step = in->step;
    double *S_forw_in = in->S_forw;
    
    int num_stages = args->num_stages;
    int newton_iter = args->newton_iter;
    double *A_mat = args->A_mat;
    double *b_vec = args->b_vec;
 
	double *rGt = workspace->rGt; 
	double *jac_out = workspace->jac_out;
    double *Jt = workspace->Jt;
    double *ode_args = workspace->ode_args;
    int *ipiv = workspace->ipiv;
	struct blasfeo_dmat *JGK = workspace->JGK;
	struct blasfeo_dvec *rG = workspace->rG;
    struct blasfeo_dvec *K = workspace->K;
    struct blasfeo_dmat *JGf = workspace->JGf;
    struct blasfeo_dmat *JKf = workspace->JKf;
    struct blasfeo_dvec *xt = workspace->xt;
    struct blasfeo_dvec *xn = workspace->xn;
    struct blasfeo_dmat *S_forw = workspace->S_forw;

    // for adjoint
    struct blasfeo_dvec *lambda = workspace->lambda;
    struct blasfeo_dvec *lambdaK = workspace->lambdaK;
    struct blasfeo_dvec *xn_traj = workspace->xn_traj;
    struct blasfeo_dvec *K_traj = workspace->K_traj;
    struct blasfeo_dmat *JG_traj = workspace->JG_traj;
    double *S_adj_in = workspace->S_adj_w;

    double *x_out = out->xn;
    double *S_forw_out = out->S_forw;
    double *S_adj_out = out->S_adj;

    acados_timer timer, timer_ad;
    double timing_ad = 0.0;

    // initialize
    blasfeo_dgese(nx*num_stages, nx*num_stages, 0.0, JGK, 0, 0);
    blasfeo_dgese(nx*num_stages, nx+nu, 0.0, JGf, 0, 0);
    blasfeo_dgese(nx*num_stages, nx+nu, 0.0, JKf, 0, 0);
    blasfeo_pack_dmat(nx, nx+nu, S_forw_in, nx, S_forw, 0, 0);

    blasfeo_dvecse(nx*num_stages, 0.0, rG, 0);
    blasfeo_dvecse(nx*num_stages, 0.0, K, 0);
    blasfeo_pack_dvec(nx, x, xn, 0);

    blasfeo_dvecse(nx*num_stages, 0.0, lambdaK, 0);

    for(kk=0;kk<nx;kk++)
        S_adj_in[kk] = in->S_adj[kk];
    for(kk=0;kk<nu;kk++)
        S_adj_in[nx+kk] = 0.0;
    blasfeo_pack_dvec(nx+nu, S_adj_in, lambda, 0);

    for (kk=0;kk<2*nx;kk++) //initialize x,xdot with zeros
        ode_args[kk] = 0.0;
    for (kk=0;kk<nu;kk++) //set controls
        ode_args[2*nx+kk] = u[kk];
    for (kk=0;kk<np;kk++)  // set parameters
        ode_args[2*nx+nu+kk] = p[kk];

    for (kk=0;kk<nx;kk++){
        rGt[kk] = 0.0;
    }
    for (kk=0;kk<nx*(2*nx+nu);kk++)
        jac_out[kk] = 0.0;
    for (kk=0;kk<nx*nx;kk++)
        Jt[kk] = 0.0;
    
        
    // start the loop
    acados_tic(&timer);
    for(ss=0; ss<num_steps; ss++){

        //  obtain Kn
        for(iter=0; iter<newton_iter; iter++){

            if (args->sens_adj){
                blasfeo_dveccp(nx, xn, 0, &xn_traj[ss], 0);
            }
           
            for(ii=0; ii<num_stages; ii++){ // ii-th row of tableau   
                    
                // take x(n); copy a strvec into a strvec
                blasfeo_dveccp(nx, xn, 0, xt, 0);

                for(jj=0; jj<num_stages; jj++){ // jj-th col of tableau
                    a = A_mat[ii+num_stages*jj];
                    if(a!=0){           // xt = xt + T_int * a[i,j]*K_j
                        a *= step; 
                        blasfeo_daxpy(nx, a, K, jj*nx, xt, 0, xt, 0); 
                    }
                }
                // put xn+sum kj into first nx elements of ode_arg               
                blasfeo_unpack_dvec(nx, xt, 0, ode_args);
                    
                // put ki into the next nx elements of ode_args
                blasfeo_unpack_dvec(nx, K, ii*nx, ode_args+nx);      

                // compute the residual of implicit ode at time t_ii, store value in rGt
                compute_impl_res(nx, nu, np, ode_args, rGt, args, mem, workspace);
                // fill in elements of rG  - store values rGt on (ii*nx)th position of rG
                blasfeo_pack_dvec(nx, rGt, rG, ii*nx);
                    
                if ( (ss==0 && iter==0 && args->jac_reuse) || (!args->jac_reuse) ){
                    // compute the jacobian of implicit ode
                    acados_tic(&timer_ad);
                    compute_impl_jac(nx, nu, np, ode_args, jac_out, args, mem, workspace);
                    timing_ad += acados_toc(&timer_ad);

                    // compute the blocks of JGK
                    for (jj=0; jj<num_stages; jj++){ //compute the block (ii,jj)th block = Jt
                        a = A_mat[ii + num_stages*jj];
                        if (a!=0){
                            a *= step;
                            for (kk=0; kk<nx*nx; kk++)
                                Jt[kk] = a* jac_out[kk];                          
                        }
                        if(jj==ii){
                            for (kk=0;kk<nx*nx;kk++)
                                Jt[kk] += jac_out[nx*nx+kk];
                        }
                        // fill in the ii-th, jj-th block of JGK
                        blasfeo_pack_dmat(nx, nx, Jt, nx, JGK, ii*nx, jj*nx);
                    } // end jj
                }
            } // end ii
                
            //DGETRF computes an LU factorization of a general M-by-N matrix A
            //using partial pivoting with row interchanges.
            if ( (ss==0 && iter==0 && args->jac_reuse) || (!args->jac_reuse) ){
                blasfeo_dgetrf_rowpivot(nx*num_stages, nx*num_stages, JGK, 0, 0, JGK, 0, 0, ipiv);
            }
        
            // permute also the r.h.s
            blasfeo_dvecpe(nx*num_stages, ipiv, rG, 0);

            // solve JGK * y = rG, JGK on the (l)eft, (l)ower-trian, (n)o-trans
            //                    (u)nit trian
            blasfeo_dtrsv_lnu(nx*num_stages, JGK, 0, 0, rG, 0, rG, 0);
            
            // solve JGK * x = rG, JGK on the (l)eft, (u)pper-trian, (n)o-trans
            //                    (n)o unit trian , and store x in rG
            blasfeo_dtrsv_unn(nx*num_stages, JGK, 0, 0, rG, 0, rG, 0);
            // scale and add a generic strmat into a generic strmat // K = K - rG, where rG is DeltaK
            blasfeo_daxpy(nx*num_stages, -1.0, rG, 0, K, 0, K, 0);
        }// end iter

        if (args->sens_adj){
            blasfeo_dveccp(nx*num_stages, K, 0, &K_traj[ss], 0);
        }
        
        // evaluate forward sensitivities
        if (args->sens_forw){

            // evaluate JGK(xn,Kn)         
            for(ii=0; ii<num_stages; ii++){  

                blasfeo_dveccp(nx, xn, 0, xt, 0);
                //compute xt = final x;
                for(jj=0; jj<num_stages; jj++){ 
                    a = A_mat[ii+num_stages*jj];
                    if(a!=0){
                        a *= step;
                        blasfeo_daxpy(nx, a, K, jj*nx, xt, 0, xt, 0);
                    }
                }               
                blasfeo_unpack_dvec(nx, xt, 0, ode_args);                                           
                blasfeo_unpack_dvec(nx, K, ii*nx, ode_args+nx);  
                
                acados_tic(&timer_ad);
                compute_impl_jac(nx, nu, np, ode_args, jac_out, args, mem, workspace);
                blasfeo_pack_dmat(nx, nx, jac_out, nx, JGf, ii*nx, 0);
                blasfeo_pack_dmat(nx, nu, jac_out+2*nx*nx, nx, JGf, ii*nx, nx);

                timing_ad += acados_toc(&timer_ad);
                
                for (jj=0;jj<num_stages;jj++){
                    a = A_mat[ii+num_stages*jj];
                    if (a!=0){
                        a *= step;
                        for (kk=0;kk<nx*nx;kk++)
                            Jt[kk] = a* jac_out[kk];                          
                    }
                    if(jj==ii){
                        for (kk=0;kk<nx*nx;kk++)
                            Jt[kk] += jac_out[nx*nx+kk];
                    }
                    blasfeo_pack_dmat(nx, nx, Jt, nx, JGK, ii*nx, jj*nx);            
                } // end jj     
            } // end ii

            // factorize JGK
            blasfeo_dgetrf_rowpivot(nx*num_stages, nx*num_stages, JGK, 0, 0, JGK, 0, 0, ipiv+(ss+1)*num_stages*nx);
            for (jj=0;jj<nx*num_stages;jj++)
                ipiv[jj] = ipiv[(ss+1)*num_stages*nx + jj];

            if (args->sens_adj){ // store the factorization and permutation
                blasfeo_dgecp(nx*num_stages, nx*num_stages, JGK, 0, 0, &JG_traj[ss], 0, 0);
            }

            // obtain JKf
            blasfeo_dgemm_nn(nx*num_stages, nx, nx, 1.0, JGf, 0, 0, S_forw, 0, 0, 0.0, JKf, 0, 0, JKf, 0, 0);
            blasfeo_dgemm_nn(nx*num_stages, nu, nx, 1.0, JGf, 0, 0, S_forw, 0, nx, 1.0, JGf, 0, nx, JKf, 0, nx);  
            blasfeo_drowpe(nx*num_stages, ipiv, JKf);
            blasfeo_dtrsm_llnu(nx*num_stages, nx+nu, 1.0, JGK, 0, 0, JKf, 0, 0, JKf, 0, 0); 
            blasfeo_dtrsm_lunn(nx*num_stages, nx+nu, 1.0, JGK, 0, 0, JKf, 0, 0, JKf, 0, 0);

            // update forward sensitivity
            for(jj=0; jj<num_stages; jj++) {
                blasfeo_dgead(nx, nx+nu, -step*b_vec[jj], JKf, jj*nx, 0, S_forw, 0, 0);
            }
        }

        // obtain x(n+1)
        for(ii=0;ii<num_stages;ii++){
            blasfeo_daxpy(nx, step*b_vec[ii], K, ii*nx, xn, 0, xn, 0);         
        }
            
    }// end int step ss

    // evaluate backwards
    if (args->sens_adj){ 
              
        for(ss=num_steps-1;ss>-1;ss--){

            blasfeo_dveccp(nx, &xn_traj[ss], 0, xt, 0);

            if (args->sens_forw){ // evalute JGf and extract factorization

                for(ii=0; ii<num_stages; ii++){  

                    for(jj=0; jj<num_stages; jj++){ 
                        a = A_mat[ii+num_stages*jj];
                        if(a!=0){
                            a *= step;
                            blasfeo_daxpy(nx, a, &K_traj[ss], jj*nx, xt, 0, xt, 0);
                        }
                    }               
                    blasfeo_unpack_dvec(nx, xt, 0, ode_args);                                           
                    blasfeo_unpack_dvec(nx, &K_traj[ss], ii*nx, ode_args+nx);  

                    acados_tic(&timer_ad);
                    compute_impl_jac(nx, nu, np, ode_args, jac_out, args, mem, workspace);
                    blasfeo_pack_dmat(nx, nx, jac_out, nx, JGf, ii*nx, 0);
                    blasfeo_pack_dmat(nx, nu, jac_out+2*nx*nx, nx, JGf, ii*nx, nx);
                    timing_ad += acados_toc(&timer_ad);
                }
            
            }else{ // extract trajectory to evaluate JGf and JGK

                for(ii=0; ii<num_stages; ii++){  

                    for(jj=0; jj<num_stages; jj++){ 
                        a = A_mat[ii+num_stages*jj];
                        if(a!=0){
                            a *= step;
                            blasfeo_daxpy(nx, a, &K_traj[ss], jj*nx, xt, 0, xt, 0);
                        }
                    }               
                    blasfeo_unpack_dvec(nx, xt, 0, ode_args);                                           
                    blasfeo_unpack_dvec(nx, &K_traj[ss], ii*nx, ode_args+nx);  

                    acados_tic(&timer_ad);
                    compute_impl_jac(nx, nu, np, ode_args, jac_out, args, mem, workspace);
                    blasfeo_pack_dmat(nx, nx, jac_out, nx, JGf, ii*nx, 0);
                    blasfeo_pack_dmat(nx, nu, jac_out+2*nx*nx, nx, JGf, ii*nx, nx);
                    timing_ad += acados_toc(&timer_ad);
                    for (jj=0;jj<num_stages;jj++){
                        a = A_mat[ii+num_stages*jj];
                        if (a!=0){
                            a *= step;
                            for (kk=0;kk<nx*nx;kk++)
                                Jt[kk] = a* jac_out[kk];                          
                        }
                        if(jj==ii){
                            for (kk=0;kk<nx*nx;kk++)
                                Jt[kk] += jac_out[nx*nx+kk];
                        }
                        blasfeo_pack_dmat(nx, nx, Jt, nx, &JG_traj[ss], ii*nx, jj*nx);            
                    } // end jj     
                } // end ii

                // factorize JGK
                blasfeo_dgetrf_rowpivot(nx*num_stages, nx*num_stages, &JG_traj[ss], 0, 0, &JG_traj[ss], 0, 0, ipiv+(ss+1)*num_stages*nx); //
            }// else if/else

            for(jj=0; jj<num_stages; jj++) {
                blasfeo_dveccpsc(nx, -step*b_vec[jj], lambda, 0, lambdaK, jj*nx);
            }

            blasfeo_dtrsv_utn(nx*num_stages, &JG_traj[ss], 0, 0, lambdaK, 0, lambdaK, 0);

            blasfeo_dtrsv_ltu(nx*num_stages, &JG_traj[ss], 0, 0, lambdaK, 0, lambdaK, 0);

            blasfeo_dvecpei(nx*num_stages, ipiv+(ss+1)*num_stages*nx, lambdaK, 0);

            blasfeo_dgemv_t(nx*num_stages, nx+nu, 1.0, JGf, 0, 0, lambdaK, 0, 1.0, lambda, 0, lambda, 0);
        }
    }

    // extract output

    blasfeo_unpack_dvec(nx, xn, 0, x_out);

    blasfeo_unpack_dmat(nx, nx+nu, S_forw, 0, 0, S_forw_out, nx);

    blasfeo_unpack_dvec(nx+nu, lambda, 0, S_adj_out);

    out->info->CPUtime = acados_toc(&timer);
    out->info->LAtime = 0.0;
    out->info->ADtime = timing_ad;

    return 0;
}
